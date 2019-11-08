#!/usr/python
#
# Copyright (c) 2019 AG Freund, Freie Universitaet Berlin (GER)
#
# This program is distributed as a free scientific tool,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
#
# For questions about running the script or for reporting bugs, please contact either:
# Christian Freund (chfreund@zedat.fu-berlin.de)
# or
# Eliot Morrison (eliot.morrison@fu-berlin.de)
#

import sys
import re
import csv
import os
import numpy as np
from scipy import stats
import math
import itertools
from pyteomics import fasta, parser, mass, achrom, electrochem, auxiliary
from scipy import stats
from os import listdir
from os.path import isfile, join
from scipy.stats import ttest_ind
import statistics

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
from matplotlib_venn import venn3, venn3_circles, venn2, venn2_circles
import matplotlib.patches as patches

import cgi
import cgitb
#cgitb.enable()
import timeit
import html

from PIL import Image, ImageDraw, ImageFont


#expdir = sys.argv[7]
#resultdir = '../nobackup/results/'+(expdir.split('../../nobackup/uploads/')[-1])
expdir = './'
resultdir = './'
net_path = '/Users/eliotmorrison/netMHCIIpan-3.2/'
iedb_path = '../IEDB/'


if not os.path.exists(resultdir):
    os.mkdir(resultdir, 0o755 )
else:
    pass

def openfile(filename, cols):
    filename = str(filename)
    filename = re.sub("b'", "", filename)
    filename = re.sub("'", "", filename)
    if 'final_renorm' not in filename:
        filename = os.path.expanduser(expdir+filename.split('/')[-1])
    else:
        filename = os.path.expanduser(resultdir+filename.split('/')[-1])
	
    raw = list(csv.reader(open(filename, 'r'), delimiter='\t'))
    headers = raw[0]
    raw = raw[1:]

    count = 0
    inds = []
    for a in headers:
        for b in cols:
            if a == b:
                inds.append(count)
        count += 1

    return raw, headers, inds

def savefile(filename, out, headers_out):
    filename = str(filename)
    filename = re.sub("b'", "", filename)
    filename = re.sub("'", "", filename)
    if headers_out != []:
        out.insert(0,headers_out)

    filename = filename.split('/')[-1]
    if 'volcano_table' in filename or 'final_renorm' in filename:
        path = resultdir
    else:
        path = expdir

    with open(os.path.expanduser(path+filename), 'w') as a_file:
        for result in out:
            result = '\t'.join(result)
            a_file.write(result + '\n')



def get_geo_mean(ratios):
    rat_sum = 0.0
    rat_prod = 1.0
    for b in ratios:
        rat_sum += b
        rat_prod *= b
    avg = float(rat_sum) / float(len(ratios))
    geo_mean = (rat_prod)**(1.0/float(len(ratios)))
    geo_sum = 0.0
    for i in ratios:
        sqr = math.log(float(i)/float(geo_mean))**2
        geo_sum += sqr
    sd_geo = math.exp(math.sqrt(geo_sum/float(len(ratios))))
    # SD(geo) is never less than 1.0. It is a factor.
    # If the geo mean is 0.91, and SD(geo) is 1.06, then the +/- SD confidence interval 
    # is 0.86 (which is 0.91 / 1.06) to 0.96 (0.91*1.06). The 95% confidence interval is 
    # 0.91 / 1.06**(1.96) = 0.81 to 0.91 * 1.06**1.96 = 1.02.
    # - 95% confidence interval
    
    minus_95_conf = geo_mean / (sd_geo**(1.96))
    plus_95_conf = geo_mean * (sd_geo**(1.96))

    est = math.log(geo_mean)
    l = math.log(minus_95_conf)
    u = math.log(plus_95_conf)
    se = (u - l) / (2*1.96)
    z = abs(est / se)
    p = math.exp(-0.717*z - 0.416*(z**2))

    sem = statistics.stdev(ratios)/math.sqrt(float(len(ratios)))

    return geo_mean, sd_geo, p, sem





def html_output(content):
    # Headers
    #print("Content-Type:text/html; charset=UTF-8")
    #print() # blank line required at end of headers

    print('<center>')
    # Content
    #print("<html><body>")
    for x in content:
        print(x + '<br>')
    #print("</body></html>")



# This gets the total intensity areas for each raw file
def add_areas_to_evidence(evidence_file):
    cols = ['Raw file', 'm/z', 'Retention time', 'Intensity']
    raw, headers, inds = openfile(evidence_file, cols)

    unique_exps = []
    headers_out = ['Run', 'Total intensity']
    out = []

    for a in raw:
        unique_exps.append(a[inds[0]])

    unique_exps = list(set(unique_exps))
    unique_exps.sort()

    for a in unique_exps:
        entry = ['']*len(headers_out)
        entry[0] = a

        intens_sum = 0.0
        for b in raw:
            if b[inds[0]] == entry[0] and b[inds[3]] != '':
                intens_sum += float(b[inds[3]])
        entry[1] = str(intens_sum)
        out.append(entry)

    out_file = 'ms_exps_total_intens_areas_for_norm.txt'
    savefile(out_file, out, headers_out)

    exp_tots = out[1:]

    headers_out = headers
    headers_out.append('Intensity')
    headers_out.append('Intensity (normalized)')

    out = []
    for a in raw:
        entry = ['']*len(headers_out)
        for x in range(0,len(a)):
            entry[x] = a[x]
        last_ind = x + 1

        intens = a[inds[3]]
        exp = a[inds[0]]

        for b in exp_tots:
            if b[0] == exp:
                intens_tot = float(b[1])

        intens_norm = ''
        if intens != '':
            intens_norm = float(intens)/intens_tot

        entry[last_ind] = str(intens)
        entry[last_ind+1] = str(intens_norm)
        out.append(entry)

    evidence_file = str(evidence_file)
    out_file = evidence_file.split('.txt')[0] + '_intens_norm.txt'
    savefile(out_file, out, headers_out)








def get_cond_peps(file_in, fasta_file):
    cols = ['Sequence', 'Proteins', 'Raw file', 'Intensity (normalized)']
    raw, headers, inds = openfile(file_in, cols)

    unique_prots_raw = []
    unique_conds = []
    unique_bio_reps = ['']
    count = 0
    prot_seqs = []
    no_fasta = []

    for a in raw:
        prots = a[inds[1]].split(';')
        for b in prots:
            if len(b) > 0:
                prot = b.split('CON__')[-1]
                if '|' in prot:
                    prot = prot.split('|')[1]
                unique_prots_raw.append(prot)
        unique_conds.append(a[inds[2]])

    unique_prots_raw = list(set(unique_prots_raw))
    unique_conds = list(set(unique_conds))
    unique_prots = []

    fasta_str = ''
    fasta_seqs = []
    for description, sequence in fasta.read(fasta_file):
        if '|' in description:
            fasta_str += description.split('|')[1] + ';'
        else:
            fasta_str += description + ';'
        fasta_seqs.append(sequence)

    count = 0
    for prot in unique_prots_raw:
        count += 1
        fasta_test = False
        if prot+';' in fasta_str:
            fasta_test = True
            #print('Fasta match found:', prot, count, '/', len(unique_prots_raw))
            unique_prots.append(prot)
        else:
            no_fasta.append([prot])

    unique_prots = list(set(unique_prots))

    headers_out = ['Protein', 'Experiment', 'Passing Peptides', 'FASTA seq', 'Intensities (normalized)']
    out = []
	
    count = 0

    small_fasta = []

    for prot in unique_prots:
        count += 1
        split1 = fasta_str.split(prot+';')
        fasta_ind = len(split1[0].split(';'))
        prot_seq = fasta_seqs[fasta_ind-1]

        small_fasta.append(['>|' + prot + '|'])
        small_fasta.append([prot_seq])

        prot_raw = []
        for a in raw:
            if prot in a[inds[1]]:
                prot_raw.append(a)

        for cond in unique_conds:
            pep_str = ''
            intens_str = ''
            for a in prot_raw:
                if cond == a[inds[2]]:
                    pep_str += a[inds[0]] + ';'
                    intens_str += a[inds[3]] + ';'

            pep_str = pep_str[:-1]
            intens_str = intens_str[:-1]

            if pep_str != '':
                entry = [prot, cond, pep_str, prot_seq, intens_str]
                #print(prot, count, '/', len(unique_prots))
                out.append(entry)
												
    file_out = file_in.split('.txt')[0] + '_passpeps.txt'
    savefile(file_out, out, headers_out)

    file_out = file_in.split('.txt')[0] + '_missing_fasta.txt'
    savefile(file_out, no_fasta, ['Missing'])

    file_out = fasta_file.split('.fasta')[0] + '_small.fasta'
    savefile(file_out, small_fasta, [])


def get_cond_peps_filt(file_in, fasta_file, var_file, bio_rep_min, tech_rep_min):
    cols = ['Sequence', 'Proteins', 'Raw file', 'Raw file', 'Intensity (normalized)']
    raw, headers, inds = openfile(file_in, cols)

    unique_peps = []
    for a in raw:
        unique_peps.append(a[inds[0]])

    unique_peps = list(set(unique_peps))

    unique_conds_1 = []
    unique_bio_reps_1 = []
    unique_tech_reps_1 = []

    var_raw, v_h, v_i = openfile(os.path.expanduser(var_file), [])
    vars = []
    var_raw.insert(0,v_h)
    for a in var_raw:
        raw_file = a[0]
        cond = a[1]
        bio_rep = a[2]
        tech_rep = a[3]

        unique_conds_1.append(cond)
        unique_bio_reps_1.append(cond+';'+bio_rep)
        unique_tech_reps_1.append(cond+';'+bio_rep+';'+tech_rep)
        vars.append(a)

    unique_conds_1 = list(set(unique_conds_1))
    unique_bio_reps_1 = list(set(unique_bio_reps_1))
    unique_tech_reps_1 = list(set(unique_tech_reps_1))

    pep_headers = ['Sequence', 'Protein']
    for cond in unique_conds_1:
        pep_headers.append('Pass in ' + cond)
    pep_pass_table = []

    unique_peps.sort()
    raw.sort()
    start_ind = 0

    count = 0
    for pep in unique_peps:
        count += 1
        pass_test = False
        pep_ind = 2
        entry = ['']*len(pep_headers)
        entry[0] = pep
        prots = []

        pep_raw = []
        for i in range(start_ind,len(raw)-1):
            a = raw[i]
            next_ent = raw[i+1]
            if i >= len(raw)-1:
                pep_raw.append(a)

            if a[inds[0]] == pep:
                pep_raw.append(a)
                if next_ent[inds[0]] != pep:
                    start_ind = i-1
                    break

        for cond in unique_conds_1:
            cond_raws = []
            for b in vars:
                if b[1] == cond:
                    cond_raws.append(b[0])
	
            bio_pass = False
            bio_reps_found = []
            for a in pep_raw:
                if a[inds[0]] == pep:
                    prots.append(a[inds[1]])
                    check = False
                    raw_file = a[inds[2]]
                    for c in cond_raws:
                        if c == raw_file:
                            check = True
                            for d in vars:
                                if d[0] == raw_file:
                                    bio_rep = d[2]
                                    tech_rep = d[3]
                    if check == True:
                        if bio_rep not in bio_reps_found:
                            bio_reps_found.append(bio_rep)

            if len(bio_reps_found) >= bio_rep_min:
                bio_nums = []
                for bio_rep_a in bio_reps_found:
                    bio_rep_raws = []
                    for c in vars:
                        if c[1] == cond and c[2] == bio_rep_a:
                            bio_rep_raws.append(c[0])

                    tech_reps_found = []
                    for a in pep_raw:
                        if a[inds[0]] == pep:
                            check = False
                            raw_file = a[inds[2]]
                            for c in bio_rep_raws:
                                if c == raw_file:
                                    check = True
                                    for d in vars:
                                        if d[0] == raw_file:
                                            bio_rep = d[2]
                                            tech_rep = d[3]
                            if check == True and bio_rep_a == bio_rep:
                                tech_reps_found.append(tech_rep)
                    tech_reps_found = list(set(tech_reps_found))
                    bio_nums.append([bio_rep_a, len(tech_reps_found)])

                passing_bios = 0
                for d in bio_nums:
                    if d[1] >= tech_rep_min:
                        passing_bios += 1

                if passing_bios >= bio_rep_min:
                    entry[pep_ind] = 'Yes'
                else:
                    entry[pep_ind] = 'No'

            else:
                entry[pep_ind] = 'No'				

            pep_ind += 1

        prots = list(set(prots))
        prot_str = ''
        for b in prots:
            prot_str += b + ';'
        entry[1] = prot_str
        #print(pep, count, '/', len(unique_peps))

        pep_pass_table.append(entry)


    pep_file_out = file_in.split('.txt')[0] + '_pep_pass_table.txt'
    savefile(pep_file_out, pep_pass_table, pep_headers)
    pep_pass_table = pep_pass_table[1:]

    unique_prots_raw = []
    unique_conds = []
    unique_bio_reps = ['']
    count = 0
    prot_seqs = []
    no_fasta = []

    for a in raw:
        prots = a[inds[1]].split(';')
        count += 1
        for b in prots:
            if len(b) > 0:
                prot = b.split('CON__')[-1]
                if '|' in prot:
                    prot = prot.split('|')[1]
                unique_prots_raw.append(prot)
        unique_conds.append(a[inds[2]])

    unique_prots_raw = list(set(unique_prots_raw))
    unique_conds = list(set(unique_conds))
    unique_prots = []

    fasta_str = ''
    fasta_seqs = []
    for description, sequence in fasta.read(fasta_file):
        if '|' in description:
            fasta_str += description.split('|')[1] + ';'
        else:
            fasta_str += description + ';'
        fasta_seqs.append(sequence)

    count = 0
    for prot in unique_prots_raw:
        count += 1
        fasta_test = False
        if prot+';' in fasta_str:
            fasta_test = True
            unique_prots.append(prot)
        else:
            no_fasta.append([prot])
    unique_prots = list(set(unique_prots))

    headers_out = ['Protein', 'Experiment', 'Passing Peptides', 'FASTA seq', 'Intensities (normalized)']
    out = []

    small_fasta = []

    count = 0
    for prot in unique_prots:
        count += 1
        split1 = fasta_str.split(prot+';')
        fasta_ind = len(split1[0].split(';'))
        prot_seq = fasta_seqs[fasta_ind-1]
        
        small_fasta.append(['>|' + prot + '|'])
        small_fasta.append([prot_seq])

        prot_raw = []
        for a in raw:
            if prot in a[inds[1]]:
                prot_raw.append(a)
        pep_raw = []
        for a in pep_pass_table:
            if prot in a[1]:
                pep_raw.append(a)

        ents = []
        zeroes = 0
        for raw_file in unique_conds:
            pep_str = ''
            intens_str = ''
            for a in prot_raw:
                if raw_file == a[inds[2]]:
                    pass_test = False
                    for b in vars:
                        if raw_file == b[0]:
                            cond = b[1]
                    for i in range(0,len(pep_headers)):
                        if 'Pass in ' + cond in pep_headers[i]:
                            pep_ind = i
                    for c in pep_raw:
                        if c[0] == a[inds[0]]:
                            if c[pep_ind] == 'Yes':	
                                pep_str += a[inds[0]] + ';'
                                intens_str += a[inds[4]] + ';'

            pep_str = pep_str[:-1]
            intens_str = intens_str[:-1]

            entry = [prot, raw_file, pep_str, prot_seq, intens_str]
            ents.append(entry)
            if intens_str == '':
                zeroes += 1
            #out.append(entry)
        if zeroes != len(ents):
            for entry in ents:
                out.append(entry)
                print(entry[-1], count, '/', len(unique_prots))
												
    file_out = file_in.split('.txt')[0] + '_passpeps_filt.txt'
    savefile(file_out, out, headers_out)

    file_out = file_in.split('.txt')[0] + '_missing_fasta.txt'
    savefile(file_out, no_fasta, ['Missing'])

    file_out = fasta_file.split('.fasta')[0] + '_small.fasta'
    savefile(file_out, small_fasta, [])


def gen_epitopes(file_in, fasta_file, min_epi_len, min_step_size, min_epi_overlap):
    cols = ['Protein', 'Passing Peptides', 'FASTA seq', 'Intensities (normalized)']
    raw, headers, inds = openfile(file_in, cols)

    headers_out = headers[:-1]
    headers_out.append('Core epitopes')
    headers_out.append('Whole epitopes')
    headers_out.append('PLAtEAU Intens. (norm)')
    headers_out.append('Peptides contributing to core epitopes')
    out = []

    jump_count = 0
    total_count = 0
						
    for j in raw:
        total_count += 1
        entry = ['']*len(headers_out)
        for i in range(0,len(j)-1):
            entry[i] = j[i]
        core_epis = []
        whole_epis = []
        founds_for_entry = []
        pep_matches = []

        prot = j[inds[0]]
        prot_seq = j[inds[2]]

        if j[inds[1]] != '':
            pep_list = j[inds[1]].split(';')
            intens_list = j[inds[3]].split(';')
            # 0 = AA, 1 = position
            seq_pos = []
            count = 1
            for a in prot_seq:
                seq_pos.append([a, count])
                count += 1

            # pep_pos: 0 = peptide, 1 = start_index
            pep_pos = []

            # 0 = AA, 1 = position
            tot_found = []

            for y in range(0,len(pep_list)):
                pep = pep_list[y]
                intens = intens_list[y]
                if intens == '':
                    intens = '0.0'
								
                before_pep = prot_seq.split(pep)[0]
                if len(prot_seq.split(pep)) > 1:
                    after_pep = prot_seq.split(pep)[1]
                else:
                    after_pep = ''

                pep_start_ind = len(before_pep)
                pep_pos.append([pep, pep_start_ind])

                for i in range(0,len(pep)):
                    tot_found.append([pep[i], pep_start_ind + i + 1, float(intens)])

            pep_pos.sort(key=lambda x: int(x[-1]))

            # this gives a total protein primary sequence w/ spec. count at each position
            seq_list = []
            for a in seq_pos:
                tot_count = 0 
                intens_count = 0.0
                for b in tot_found:
                    if b[0] == a[0] and b[1] == a[1]:
                        tot_count += 1
                        intens_count += b[2]

                #seq_list.append([a[0], a[1], tot_count])
                seq_list.append([a[0], a[1], intens_count])

            # Trim the seq_list to 5 residues before first intens
            # and 5 residues after last intens, to minimize time
            # spent going over 0-intensity residues

            epis = []
            start_pos = []
            end_pos = []
            for i in range(1,len(seq_list)-1):
                prev_num = float(seq_list[i-1][2])
                current_num = float(seq_list[i][2])
                next_num = float(seq_list[i+1][2])

                #if current_num != 0 or next_num != 0:
                    #print(i, prev_num, current_num, next_num)

                if current_num == 0.0 and next_num != 0.0:
                    #start_pos.append(i-3) # ******
                    start_pos.append(i) # ******
                    #print('Start pos:', i+1)
                if i == 1 and current_num != 0.0:
                    start_pos.append(i-1) 
                    #print('Start pos:', i+1)

                if current_num != 0.0 and next_num == 0.0:
                    end_pos.append(i+2)
                    #print('End pos:', i)
                    #print('')
                if i == len(seq_list)-2 and current_num != 0.0:
                    end_pos.append(i+2)


            new_start_pos = []
            new_end_pos = []

            for i in range(0,len(start_pos)):
                epi = seq_list[start_pos[i]:end_pos[i]]

                found_str = ''
                founds = []
                for x in epi:
                    if float(x[2]) not in founds:
                        if float(x[2]) != 0.0:
                            founds.append(float(x[2]))
                    found_str += str(x[2]) + '--'
                found_str = found_str[:-2]

                # a normal epitope should go low > high > low
                # two local maxima will go low > high > low > high > low

                # Create a string: 0 = less than or equal, 1 = greater than
                found_strs = []
                founds.sort()
                for y in founds:
                    bin_str = ''
                    for z in epi:
                        if float(z[2]) <= y:
                            bin_str += '0'
                        else:
                            bin_str += '1'
                    #print(y, bin_str)
                    split_1 = bin_str.split('0')
                    split_2 = []
                    for z in split_1:
                        if z != '':
                            split_2.append(z)
                    #print(y, split_1, split_2)							
                    if len(split_2) > 1:
                        if len(split_2[0]) >= min_epi_len and len(split_2[-1]) >= min_epi_len:
                            found_strs.append([y, len(split_2[0]), len(split_2[-1]), bin_str])
                            #print(y, bin_str)
                            #print(found_str)
                            #print(found_strs[-1])
                if len(found_strs) > 0:
                    #print(found_strs)
                    new_min = found_strs[0][0]
                    bin_str = found_strs[0][-1]
                    zero_list = []
                    for k in bin_str.split('1'):
                        if k != '':
                            zero_list.append(k)
                    gap_size = len(zero_list[1])
                    gap_pos_bin = len(zero_list[0]) + found_strs[0][1]
                    gap_pos = epi[gap_pos_bin][1]
                    #print('Gap size:', gap_size, 'starting:', gap_pos_bin, gap_pos)

                    # gap_pos = first '0' of gap
                    # if gap size = 1, add new end = gap_pos - 1, new start = gap_pos
                    # if gap size > 1:
                    #				new end = len(first 1s) + gap_size
                    #				new start = len(first 1s)
                    new_starts = []
                    new_ends = []

                    if gap_size == 1:
                        new_ends.append(gap_pos-1)
                        new_end_pos.append(gap_pos-1)
                        new_starts.append(gap_pos)
                        new_start_pos.append(gap_pos)
												
                    elif gap_size > 1:
                        new_ends.append(gap_pos+gap_size)
                        new_end_pos.append(gap_pos+gap_size)
                        new_starts.append(gap_pos)
                        new_start_pos.append(gap_pos)
												

            if len(new_start_pos) > 0:
                for f in range(0,len(new_start_pos)):
                    start_pos.append(new_start_pos[f])	
                    end_pos.append(new_end_pos[f])

            start_pos.sort()
            end_pos.sort()

            # Look for overlapping epitopes
            # A standard nested peptide distribution should change by 1-2 residues
            # A large jump (> 5 residues) likely means a different epitope distribution
						
            # pep_pos: 0 = peptide, 1 = start_index
            pep_pos_unique = []
            for x in pep_pos:
                if x not in pep_pos_unique:
                    pep_pos_unique.append(x)
										
            new_start_pos = []
            new_end_pos = []
            jump_test = False

            for i in range(0,len(start_pos)):
                epi = seq_list[start_pos[i]:end_pos[i]]
                start_ind = start_pos[i]
                epi_str = ''
                for y in epi:
                    epi_str += str(y[0])

                match_pep_pos = []
                for x in pep_pos_unique:
                    #if epi_str in x[0]:
                    if x[0] in epi_str:
                        match_pep_pos.append(x)
												
                found_str = ''
                founds = []
                for x in epi:
                    if float(x[2]) not in founds:
                        if float(x[2]) != 0.0:
                            founds.append(float(x[2]))
                    found_str += str(x[2]) + '--'
                found_str = found_str[:-2]
                #jump_test = False

                for i in range(0,len(match_pep_pos)-1):
                    this_pep = match_pep_pos[i]
                    next_pep = match_pep_pos[i+1]
                    step_size = next_pep[1] - this_pep[1]
                    #print('&&', this_pep, next_pep, step_size, len(this_pep[0]), step_size + min_epi_overlap) 
                    if next_pep[1] - this_pep[1] > min_step_size and len(this_pep[0]) <= step_size + min_epi_overlap:
                        #print('')
                        #print('Jump:', this_pep, next_pep)
                        jump_test = True
                        if next_pep[1]-1 not in start_pos:
                            new_start_pos.append(next_pep[1]) ####
                        if this_pep[1]+len(this_pep[0])-1 not in end_pos:
                            new_end_pos.append(this_pep[1]+len(this_pep[0]))	####
                        #print('Step size:', step_size, 'pep length:', len(this_pep[0]))
                        #print('Start added:', new_start_pos[-1])
                        #print('End added:', new_end_pos[-1])

                        # pep_list_1 will be peptides with start index <= this_pep[1]
                        # pep_list_2 will be peptides with start index > this_pep[1]
                        pep_list_1_ind_end = this_pep[1]
                        pep_list_1_ind_start = start_ind
                        pep_list_2_ind_start = next_pep[1]
                        pep_list_2_ind_end = match_pep_pos[-1][1] + len(match_pep_pos[-1][0])
                        #pep_list_2_ind = next_pep[1]
                        pep_list_1 = []
                        pep_list_2 = []

                if jump_test == True:
                    for j in match_pep_pos:
                        if pep_list_1_ind_start <= j[1] <= pep_list_1_ind_end:
                            pep_list_1.append(j[0])
                        elif pep_list_2_ind_start <= j[1] <= pep_list_2_ind_end:
                            pep_list_2.append(j[0])


            if jump_test == True:
                jump_count += 1
                #print('Jumps:', jump_count, '/', total_count, ' (' + str(float(jump_count)/float(total_count)*100.0)+'%)')

            if len(new_start_pos) > 0:
                #print('**', new_start_pos)
                #print(new_end_pos)
                for f in range(0,len(new_start_pos)):
                    start_pos.append(new_start_pos[f])	
                    end_pos.append(new_end_pos[f])

            start_pos.sort()
            end_pos.sort()

            for i in range(0,len(start_pos)):
                epi = seq_list[start_pos[i]:end_pos[i]]
                core_start = start_pos[i]
                core_end = end_pos[i]

                founds = []
                for x in epi:
                    if float(x[2]) not in founds:
                        founds.append(float(x[2]))
                
                whole_epi = ''
                for l in epi:
                    whole_epi += l[0]


                list_1_jump_test = False
                list_2_jump_test = False
                if jump_test == True:
                    if pep_list_1_ind_start <= start_pos[i] <= pep_list_1_ind_end:
                        list_1_jump_test = True
                    if pep_list_2_ind_start <= start_pos[i] <= pep_list_2_ind_end:
                        list_2_jump_test = True

                pep_match = []
                intens_match = []
                if jump_test == False:
                    for y in range(0,len(pep_list)):
                        pep = pep_list[y]
                        intens = intens_list[y]
                        if pep in whole_epi:
                            pep_match.append(pep)
                            intens_match.append(intens)
                else:
                    #print(start_pos[i], pep_list_1_ind)
                    #print('List 1', pep_list_1_ind_start, start_pos[i], pep_list_1_ind_end)
                    #print('List 2', pep_list_2_ind_start, start_pos[i], pep_list_2_ind_end)
                    if list_1_jump_test == True:
                        for y in range(0,len(pep_list)):
                            pep = pep_list[y]
                            intens = intens_list[y]
                            if pep in whole_epi and pep in pep_list_1:
                                pep_match.append(pep)
                                intens_match.append(intens)
                    elif list_2_jump_test == True:
                        for y in range(0,len(pep_list)):
                            pep = pep_list[y]
                            intens = intens_list[y]
                            if pep in whole_epi and pep in pep_list_2:
                                pep_match.append(pep)
                                intens_match.append(intens)											

                    else:
                        for y in range(0,len(pep_list)):
                            pep = pep_list[y]
                            intens = intens_list[y]
                            if pep in whole_epi:
                                pep_match.append(pep)
                                intens_match.append(intens)

                if list_1_jump_test == True or list_2_jump_test == True:
                    tot_found_1 = []
                    for y in range(0,len(pep_match)):
                        pep = pep_match[y]
                        intens = intens_match[y]
                        if intens == '':
                            intens = 0.0
												
                        before_pep = prot_seq.split(pep)[0]
                        if len(prot_seq.split(pep)) > 1:
                            after_pep = prot_seq.split(pep)[1]
                        else:
                            after_pep = ''

                        pep_start_ind = len(before_pep)
                        #pep_pos.append([pep, pep_start_ind])

                        for n in range(0,len(pep)):
                            tot_found_1.append([pep[n], pep_start_ind + n + 1, intens])

                    # this gives a total protein primary sequence w/ spec. count at each position
                    seq_list_1 = []
                    for f in seq_pos:
                        tot_count = 0
                        intens_count = 0.0
                        for d in tot_found_1:
                            if d[0] == f[0] and d[1] == f[1]:
                                tot_count += 1
                                if d[2] != '':
                                    intens_count += float(d[2])
                        #seq_list_1.append([f[0], f[1], tot_count])
                        seq_list_1.append([f[0], f[1], intens_count])
                    epi_1 = seq_list_1[start_pos[i]:end_pos[i]]
                    founds_1 = []
                    for x in epi_1:
                        if float(x[2]) not in founds_1:
                            if float(x[2]) != 0.0:
                                founds_1.append(float(x[2]))

                pep_match_str = ''
                for g in pep_match:
                    pep_match_str += g + ','
                pep_match_str = pep_match_str[:-1]
                core_str = ''

                if list_1_jump_test == True or list_2_jump_test == True:
                    found_str = ''
                    for x in epi_1:
                        found_str += str(x[2]) + '--'
                    found_str = found_str[:-2]
                    founds_for_entry.append(found_str)

                    founds_1.sort(reverse=True)
                    passing_core = False
                    for num in founds_1:
                        core = []
                        for k in epi_1:
                            if float(k[2]) >= num:
                                core.append(k)

                        if len(core) >= min_epi_len and passing_core == False:
                            passing_core = True
                            core_str = ''

                            core_start = core[0][1]
                            core_end = core[-1][1]
                            for l in core:
                                core_str += l[0]

                    core_str = ''
                    core = seq_list[core_start-1:core_end]
                    for l in core:
                        core_str += l[0]

                    # epitope, start pos, end pos
                    if list_1_jump_test == True:
                        core_str = core_str + '*'
                    if list_2_jump_test == True:
                        core_str = '*' + core_str
                    core_epis.append([core_str, core_start, core_end])
                    if list_1_jump_test == True:
                        whole_epi = whole_epi + '*'
                    if list_2_jump_test == True:
                        whole_epi = '*' + whole_epi
                    whole_epis.append([whole_epi, start_pos[i], end_pos[i]])
                    pep_matches.append(pep_match_str)

                else:
                    found_str = ''
                    for x in epi:
                        found_str += str(x[2]) + '--'
                    found_str = found_str[:-2]
                    founds_for_entry.append(found_str)

                    founds.sort(reverse=True)
                    passing_core = False
                    for num in founds:
                        core = []
                        for k in epi:
                            if float(k[2]) >= num:
                                core.append(k)

                        if len(core) >= min_epi_len and passing_core == False:
                            passing_core = True
                            core_str = ''

                            core_start = core[0][1]
                            core_end = core[-1][1]
                            for l in core:
                                core_str += l[0]

                    core_str = ''
                    core = seq_list[core_start-1:core_end]
                    for l in core:
                        core_str += l[0]

                    # epitope, start pos, end pos
                    if list_1_jump_test == True:
                        core_str = core_str + '*'
                    if list_2_jump_test == True:
                        core_str = '*' + core_str
                    core_epis.append([core_str, core_start, core_end])
                    if list_1_jump_test == True:
                        whole_epi = whole_epi + '*'
                    if list_2_jump_test == True:
                        whole_epi = '*' + whole_epi
                    whole_epis.append([whole_epi, start_pos[i], end_pos[i]])
                    pep_matches.append(pep_match_str)
								
        core_str = ''
        for t in core_epis:
            epi_temp = str(t[1]) + '-' + t[0] + '-' + str(t[2])
            core_str += epi_temp + ';'
        core_str = core_str[:-1]
								
        whole_str = ''
        for t in whole_epis:
            epi_temp = str(t[1]) + '-' + t[0] + '-' + str(t[2])
            whole_str += epi_temp + ';'
        whole_str = whole_str[:-1]

        pep_match_str = ''
        for o in pep_matches:
            pep_match_str += o + ';'
        pep_match_str = pep_match_str[:-1]

        founds_for_entry_str = ''
        for o in founds_for_entry:
            founds_for_entry_str += o + ';'
        founds_for_entry_str = founds_for_entry_str [:-1]
				
        entry[-4] = core_str
        entry[-3] = whole_str
        entry[-2] = founds_for_entry_str
        entry[-1] = pep_match_str
        out.append(entry)
        #print(entry[0], total_count, '/', len(raw))

    file_out = file_in.split('.txt')[0] + '_epitopes.txt'
    savefile(file_out, out, headers_out)





def comb_epis(evidence_file, epitope_file, exp_name):
    cols = ['Protein', 'Experiment', 'Core epitopes', 'Whole epitopes',
            'PLAtEAU Intens. (norm)', 'Peptides contributing to core epitopes'
            ]
    epis, epi_headers, epi_inds = openfile(epitope_file, cols)

    cols = ['Sequence', 'Proteins', 'Raw file', 'Intensity (normalized)']
    raw, headers, inds = openfile(evidence_file, cols)


    # find all unique wholes
    unique_wholes = []
    unique_cores = []
    for a in epis:
        wholes = a[epi_inds[3]].split(';')
        cores = a[epi_inds[2]].split(';')
        if wholes[0] != '' and cores[0] != '':
            for i in wholes:
                whole = i.split('-')[1]
                if whole != '':
                    unique_wholes.append(whole)
            for i in cores:
                core = i.split('-')[1]
                if '--' not in i:
                    unique_cores.append(core)
    unique_wholes = list(set(unique_wholes))
    unique_cores = list(set(unique_cores))
    #print('Unique wholes:', len(unique_wholes))
    #print('Unique cores:', len(unique_cores))

    # 0 = unique core
    # 1 = longest whole match

    unique_longest_wholes = []

    for a in unique_cores:
        whole_matches = []
        whole_lens = []

        for b in unique_wholes:
            if '*' not in a and '*' not in b:
                if a in b:
                    whole_matches.append(b)
                    whole_lens.append(len(b))
            elif '*' in a:
                if a[-1] == '*' and b[-1] == '*':
                    if a[:-1] in b[:-1]:
                        whole_matches.append(b)
                        whole_lens.append(len(b))
                elif a[0] == '*' and b[0] == '*':
                    if a[1:] in b[1:]:
                        whole_matches.append(b)
                        whole_lens.append(len(b))								

        if len(whole_lens) > 0:
            longest_len = max(whole_lens)
            for i in range(0,len(whole_matches)):
                if whole_lens[i] == longest_len:
                    longest_whole = whole_matches[i]
                    if longest_whole not in unique_longest_wholes:
                        unique_longest_wholes.append(longest_whole)

    unique_longest_wholes = list(set(unique_longest_wholes))
    #print('Unique longest wholes:', len(unique_longest_wholes))

    # 0 = core, 1 = whole
    unique_longest_cores = []
    seen = []
    for a in unique_longest_wholes:
        core_matches = []
        core_lens = []
        for b in unique_cores:
            if '*' not in a and '*' not in b:
                if b in a:
                    core_matches.append(b)
                    core_lens.append(len(b))
            elif '*' in a:
                if a[-1] == '*' and b[-1] == '*':
                    if b[:-1] in a[:-1]:
                        core_matches.append(b)
                        core_lens.append(len(b))
                elif a[0] == '*' and b[0] == '*':
                    if b[1:] in a[1:]:
                        core_matches.append(b)
                        core_lens.append(len(b))

        if len(core_lens) > 0:
            longest_len = max(core_lens)
            longest_cores = []
            for i in range(0,len(core_matches)):
                if core_lens[i] == longest_len:
                    longest_cores.append(core_matches[i])

        if len(core_lens) > 0:
            longest_len = max(core_lens)
            longest_cores = []
            finals = []
            for i in range(0,len(core_matches)):
                if core_lens[i] == longest_len:
                    longest_cores.append(core_matches[i])
            if len(longest_cores) == 1:
                longest_core = longest_cores[0]
                finals.append(longest_core)
            else:
                core1 = longest_cores[0]
                core2 = longest_cores[1]
                if '*' in core1 and len(core1) > 1:
                    if core1[0] == '*':
                        core1 = core1[1:]
                    if core1[-1] == '*':
                        core1 = core1[:-1]
                if '*' in core2 and len(core2) > 1:
                    if core2[0] == '*':
                        core2 = core2[1:]
                    if core2[-1] == '*':
                        core2 = core2[:-1]

                n_dist_1 = a.split(core1)[0]
                n_dist_2 = a.split(core2)[0]
                #print(a, core1, n_dist_1)
                #print(a, core2, n_dist_2)

                if n_dist_1 > n_dist_2:
                    first_core = core2
                    second_core = core1
                elif n_dist_1 < n_dist_2:
                    first_core = core1
                    second_core = core2

                new_core = first_core + second_core[-1]
                if first_core in new_core and second_core in new_core:
                    longest_core = new_core
                    if longest_cores[0][0] == '*' and longest_cores[1][0] == '*':
                        longest_core = '*' + longest_core
                    finals.append(longest_core)
                else:
                    # two different cores of same length
                    for x in longest_cores:
                        if x not in seen:
                            seen.append(x)
                            finals.append(x)
            for longest_core in finals:
                if longest_core not in seen:
                    seen.append(longest_core)
                    unique_longest_cores.append([longest_core, a])

    #print('Unique longest cores:', len(unique_longest_cores))

    unique_conds = []
    for a in raw:
        unique_conds.append(a[inds[2]])
    unique_conds = list(set(unique_conds))
    unique_conds.sort()

    headers_out = ['Core Epitope', 'Proteins', 'Core Epitope Length', 'Whole Epitope']
    for name in unique_conds:
        headers_out.append('% Rel. Intens. in ' + name)

    out = []

    count = 0

    for a in unique_longest_cores:
        count += 1
        entry = ['']*len(headers_out)
        core = a[0]
        core_test = core
        if core[0] == '*':
            core_test = core[1:]
        if core[-1] == '*':
            core_test = core[:-1]
        whole = a[1]
        start_ind = 4

        entry[0] = core
        entry[2] = str(len(core))
        entry[3] = whole

        prots = []
        exp_intens = []

        #start1 = timeit.default_timer()
        for b in epis:
            if core_test in b[epi_inds[3]]:
                cores = b[epi_inds[2]].split(';')
                wholes = b[epi_inds[3]].split(';')
                intens = b[epi_inds[4]].split(';')

                for i in range(0,len(cores)):
                    if wholes[0] != '':
                        whole_test = wholes[i].split('-')[1]
                    else:
                        whole_test = []

                    if core_test in whole_test:
                        prots.append(b[epi_inds[0]])
                        max_intens = 0.0
                        intens_vals = []
                        for c in intens[i].split('--'):
                            intens_vals.append(float(c))
                        if len(intens_vals) > 0:
                            max_intens = max(intens_vals)
                        if max_intens != 0.0:
                            exp_intens.append([b[epi_inds[1]], max_intens])

        #stop1 = timeit.default_timer()
        #print('time: ', stop1 - start1)  

        prots = list(set(prots))
        prots.sort()
        prot_str = ''
        for b in prots:
            prot_str += b + ';'
        prot_str = prot_str[:-1]
        entry[1] = prot_str

        count_total = 0
        zeroes = 0
        all_inds = []

        for cond in unique_conds:
            count_total += 1
            intensity_sum = 0.0
            for c in exp_intens:
                if c[0] == cond:
                    intensity_sum += c[1]
            if intensity_sum == 0.0:
                zeroes += 1

            entry[start_ind] = str(intensity_sum)
            start_ind += 1

        #print(entry)
        if count_total != zeroes:
            out.append(entry)
            #print(count, '/', len(unique_longest_cores))


    file_out = str(exp_name) + '_core_epitopes_final.txt'
    savefile(file_out, out, headers_out)

				


def overlap_test_1(pep1, pep2, whole_seq):
    overlap_test = True

    split_1 = whole_seq.split(pep1)
    for x in split_1:
        if pep2 in x:
            overlap_test = False

    return overlap_test


def get_longest(group_in, whole_seq):
    seq_pos = []
    for i in range(0,len(whole_seq)):
        seq_pos.append([i, whole_seq[i]])

    core_pos = []
    for a in group_in:
        test = re.sub('\*', '', a)
        start_pos = len(whole_seq.split(test)[0])
        end_pos = start_pos + len(test)
        for i in range(start_pos, end_pos):
            core_pos.append(i)

    coverage = []
    for i in seq_pos:
        pos = i[0]
        count = 0
        for x in core_pos:
            if x == pos:
                count += 1
        coverage.append([i[0], i[1], count])
    
    start_pos = ''
    end_pos = ''

    if coverage[0][2] > 0:
        start_pos = 0
    if coverage[-1][2] > 0:
        end_pos = len(coverage) + 1
    if coverage[0][2] == 0 and coverage[1][2] > 0:
        start_pos = 1

    for i in range(1,len(coverage)-1):
        prev_pos = coverage[i-1]
        this_pos = coverage[i]
        next_pos = coverage[i+1]

        if this_pos[2] == 0 and next_pos[2] > 0:
            start_pos = i+1

        if this_pos[2] >0 and next_pos[2] == 0:
            end_pos = i+1

    core_seq = whole_seq[start_pos:end_pos]

    return core_seq






def comb_epis_2(evidence_file, epitope_file, exp_name, min_epi_len, fasta_file, epi_file):
    cols = ['Protein', 'Experiment', 'Core epitopes', 'Whole epitopes',
        'PLAtEAU Intens. (norm)', 'Peptides contributing to core epitopes'
        ]
    epis, epi_headers, epi_inds = openfile(epitope_file, cols)
    fasta_ind = 3

    cols = ['Sequence', 'Proteins', 'Experiment', 'Intensity (normalized)']
    raw, headers, inds = openfile(evidence_file, cols)

    cols = ['Protein', 'Experiment', 'Passing Peptides', 'Intensities (normalized)']
    peps, p_h, p_i = openfile(epi_file, cols)
    peps.sort(key=lambda x: x[3])

    unique_wholes = []
    unique_cores = []
    unique_fastas = []
    for a in epis:
        prots = a[epi_inds[0]]
        fasta_seq = a[fasta_ind]
        wholes = a[epi_inds[3]].split(';')
        if wholes[0] != '':
            for i in wholes:
                whole = i.split('-')[1]
                whole_test = re.sub('\*', '', whole)
                if whole_test != '' and len(whole_test) >= min_epi_len:
                    unique_wholes.append(fasta_seq + '--' + whole + '--' + prots)
                    unique_fastas.append(fasta_seq)

        cores = a[epi_inds[2]].split(';')
        if cores[0] != '':
            for i in cores:
                core = i.split('-')[1]
                if len(core) > 0:
                    if core[0] == '*':
                        core = core[2:]
                core_test = re.sub('\*', '', core)
                if core_test != '' and len(core_test) >= min_epi_len:
                    unique_cores.append(fasta_seq + '--' + core + '--' + prots)
                    unique_fastas.append(fasta_seq)

    unique_wholes = list(set(unique_wholes))
    unique_cores = list(set(unique_cores))
    unique_fastas = list(set(unique_fastas))
    unique_wholes.sort()
    unique_cores.sort()
    unique_fastas.sort()
    unique_wholes.append('--')

    unique_longest_wholes = []
    count = 0
    fasta_ind = 0
    for fasta_seq in unique_fastas:
        count += 1
        whole_matches = []
        starts = []
        ends = []
        prots = ''
        for i in range(fasta_ind,len(unique_wholes)-1):
            b = unique_wholes[i]
            next_b = unique_wholes[i+1]
            b_split = b.split('--')
            if b_split[0] == fasta_seq:
                whole_test = re.sub('\*','',b_split[1])
                start_pos = len(fasta_seq.split(whole_test)[0])
                end_pos = start_pos + len(b_split[1])
                starts.append(start_pos)
                ends.append(end_pos)
                whole_matches.append(b_split[1])
                prots = b_split[2]

                if next_b.split('--')[0] != fasta_seq:
                    fasta_ind = i-1
                    break
        whole_matches = list(set(whole_matches))

        if len(whole_matches) > 0:
            bin_seq = ''
            for i in range(0,len(fasta_seq)):
                res = 0
                for j in range(0,len(starts)):
                    if starts[j] <= i <= ends[j]:
                        res = 1
                bin_seq += str(res)

            longest_wholes = []
            ind_pos = 0
            for k in bin_seq.split('0'):
                if len(k) == 0:
                    ind_pos += 1
                else:
                    longest_wholes.append(fasta_seq[ind_pos:ind_pos+len(k)-1])
                    #print(longest_wholes[-1], whole_matches)
                    ind_pos += len(k)
                    
            longest_wholes = list(set(longest_wholes))
            print(longest_wholes, count, '/', len(unique_fastas))
            for new_whole in longest_wholes:
                unique_longest_wholes.append(fasta_seq + '--' + new_whole + '--' + prots)
                #if 'GSDQSEN' in new_whole:
                #    print(whole_matches)
                #    print(new_whole)
            
    #count = 0
    #fasta_ind = 0
    #for c in unique_cores:
    #    count += 1
    #    core = c.split('--')[1]
    #    fasta_seq = c.split('--')[0]
    #    prots = c.split('--')[2]
    #    core_test = re.sub('\*', '', core)
    #    whole_matches = []
    #    for i in range(fasta_ind,len(unique_wholes)-1):
    #        b = unique_wholes[i]
    #        next_b = unique_wholes[i+1]
    #        if b.split('--')[0] == fasta_seq:
    #            whole_test = b.split('--')[1]
    #             if i == len(unique_wholes)-1:
    #                if next_b.split('--')[0] == fasta_seq:
    #                    whole_matches.append(whole_test)
    #            if core_test in whole_test:
    #                whole_matches.append(whole_test)
    #            if next_b.split('--')[0] != fasta_seq:
    #                fasta_ind = i-1
    #                break
    #    whole_matches = list(set(whole_matches))
    #    if len(whole_matches) > 0:
    #        longest_whole = get_longest(whole_matches, fasta_seq)
    #        unique_longest_wholes.append(fasta_seq + '--' + longest_whole + '--' + prots)
    #        print('longest whole:', longest_whole, count, '/', len(unique_cores))       
    
    unique_longest_wholes = list(set(unique_longest_wholes))
    unique_longest_wholes.sort()

    test_seq = 'KKAAGAG'
    unique_longest_cores = []
    core_whole_pairs = []
    count = 0
    fasta_ind = 0
    for b in unique_longest_wholes:
        whole_seq = b.split('--')[1]
        fasta_seq = b.split('--')[0]
        prots = b.split('--')[2]
        count += 1
        #print(count, '/', len(unique_longest_wholes))
        whole_test = re.sub('\*', '', whole_seq)
        matching_cores = []

        for i in range(fasta_ind,len(unique_cores)-1):
            b = unique_cores[i]
            next_b = unique_cores[i+1]
            if b.split('--')[0] == fasta_seq:
                core_test = b.split('--')[1]
                core_test = re.sub('\*','',core_test)
                if i == len(unique_cores)-1 and next_b.split('--')[0] == fasta_seq:
                    next_test = next_b.split('--')[1]
                    if re.sub('\*','',next_test) in whole_test:
                        matching_cores.append(next_test)
                if core_test in whole_test:
                    matching_cores.append(b.split('--')[1])
                if next_b.split('--')[0] != fasta_seq:
                    fasta_ind = i-1
                    break

        matching_cores = list(set(matching_cores))

        core_pass = False
        if len(matching_cores) == 1:
            core_pass = True
            core_whole_pairs.append([matching_cores[0], whole_seq, fasta_seq, prots])
            #print(matching_cores[0])
        else:
            if len(matching_cores) == 2:
                #print(matching_cores[0], matching_cores[1])
                if '*' in matching_cores[0] and '*' in matching_cores[1]:
                    core_pass = True
                    core_whole_pairs.append([matching_cores[0], whole_seq, fasta_seq, prots])
                    core_whole_pairs.append([matching_cores[1], whole_seq, fasta_seq, prots])
                else:
                    # There are 2 cores in the same whole, but it is not a jump
                    # Check if they overlap; if so, choose the largest
                    # If not, keep them separate

                    overlap_test = overlap_test_1(matching_cores[0], matching_cores[1], whole_seq)
                    if overlap_test == False:
                        core_pass = True
                        core_whole_pairs.append([matching_cores[0], whole_seq, fasta_seq, prots])
                        core_whole_pairs.append([matching_cores[1], whole_seq, fasta_seq, prots])
                    else:
                        core_pass = True
                        longest_core = get_longest(matching_cores, whole_seq)
                        core_whole_pairs.append([longest_core, whole_seq, fasta_seq, prots])

            else:
                #print('More than 2')
                #print('')
                #print(matching_cores)
                combs = []
                for i in range(0,len(matching_cores)-1):
                    for j in range(1,len(matching_cores)):
                        test_1 = re.sub('\*', '', matching_cores[i])
                        test_2 = re.sub('\*', '', matching_cores[j])
                        overlap_test = overlap_test_1(test_1, test_2, whole_seq)
                        combs.append([matching_cores[i], matching_cores[j], overlap_test])
                        #print(combs[-1])

                group_1 = []
                group_2 = []
                group_3 = []
                group_4 = []

                for core_test in matching_cores:
                    if len(group_1) == 0:
                        group_1.append(core_test)
                        for d in combs:
                            right_star = False
                            left_star = False
                            for e in group_1:
                                if e[0] == '*':
                                    left_star = True
                                elif e[-1] == '*':
                                    right_star = True
                            if d[2] == True:
                                if d[0] == core_test:
                                    if right_star == True and d[1][0] != '*':
                                        if d[1] not in group_1:
                                            group_1.append(d[1])
                                    elif left_star == True and d[1][-1] != '*':
                                        if d[1] not in group_1:
                                            group_1.append(d[1])
                                    elif right_star == False and left_star == False:
                                        if d[1] not in group_1:
                                            group_1.append(d[1])
                                if d[1] == core_test:
                                    if right_star == True and d[0][0] != '*':
                                        if d[0] not in group_1:
                                            group_1.append(d[0])
                                    elif left_star == True and d[0][-1] != '*':
                                        if d[0] not in group_1:
                                            group_1.append(d[0])
                                    elif right_star == False and left_star == False:
                                        if d[0] not in group_1:
                                            group_1.append(d[0])
                #if test_seq in whole_test:
                    #print('Group 1:', group_1)

                for core_test in matching_cores:
                    if core_test not in group_1 and len(group_2) == 0:
                        group_2.append(core_test)
                        for d in combs:
                            right_star = False
                            left_star = False
                            for e in group_2:
                                if e[0] == '*':
                                    left_star = True
                                elif e[-1] == '*':
                                    right_star = True
                            if d[2] == True:
                                if d[0] == core_test:
                                    if right_star == True and d[1][0] != '*':
                                        if d[1] not in group_1 and d[1] not in group_2:
                                            group_2.append(d[1])
                                    elif left_star == True and d[1][-1] != '*':
                                        if d[1] not in group_1 and d[1] not in group_2:
                                            group_2.append(d[1])
                                    elif right_star == False and left_star == False:
                                        if d[1] not in group_1 and d[1] not in group_2:
                                            group_2.append(d[1])
                                if d[1] == core_test:
                                    if right_star == True and d[0][0] != '*':
                                        if d[0] not in group_1 and d[0] not in group_2:
                                            group_2.append(d[0])
                                    elif left_star == True and d[0][-1] != '*':
                                        if d[0] not in group_1 and d[0] not in group_2:
                                            group_2.append(d[0])
                                    elif right_star == False and left_star == False:
                                        if d[0] not in group_1 and d[0] not in group_2:
                                            group_2.append(d[0])
                #if test_seq in whole_test:
                    #print('Group 2:', group_2)

                for core_test in matching_cores:
                    if core_test not in group_1 and core_test not in group_2 and len(group_3) == 0:
                        group_3.append(core_test)
                        for d in combs:
                            right_star = False
                            left_star = False
                            for e in group_3:
                                if e[0] == '*':
                                    left_star = True
                                elif e[-1] == '*':
                                    right_star = True
                            if d[2] == True:
                                if d[0] == core_test:
                                    if right_star == True and d[1][0] != '*':
                                        if d[1] not in group_1 and d[1] not in group_2 and d[1] not in group_3:
                                            group_3.append(d[1])
                                    elif left_star == True and d[1][-1] != '*':
                                        if d[1] not in group_1 and d[1] not in group_2 and d[1] not in group_3:
                                            group_3.append(d[1])
                                    elif right_star == False and left_star == False:
                                        if d[1] not in group_1 and d[1] not in group_2 and d[1] not in group_3:
                                            group_3.append(d[1])
                                if d[1] == core_test:
                                    if right_star == True and d[0][0] != '*':
                                        if d[0] not in group_1 and d[0] not in group_2 and d[0] not in group_3:
                                            group_3.append(d[0])
                                    elif left_star == True and d[0][-1] != '*':
                                        if d[0] not in group_1 and d[0] not in group_2 and d[0] not in group_3:
                                            group_3.append(d[0])
                                    elif right_star == False and left_star == False:
                                        if d[0] not in group_1 and d[0] not in group_2 and d[0] not in group_3:
                                            group_3.append(d[0])
                #if test_seq in whole_test:
                    #print('Group 3:', group_3)

                for core_test in matching_cores:
                    if core_test not in group_1 and core_test not in group_2 and core_test not in group_4 and len(group_4) == 0:
                        group_4.append(core_test)
                        for d in combs:
                            right_star = False
                            left_star = False
                            for e in group_4:
                                if e[0] == '*':
                                    left_star = True
                                elif e[-1] == '*':
                                    right_star = True
                            if d[2] == True:
                                if d[0] == core_test:
                                    if right_star == True and d[1][0] != '*':
                                        if d[1] not in group_1 and d[1] not in group_2 and d[1] not in group_3 and d[1] not in group_4:
                                            group_4.append(d[1])
                                    elif left_star == True and d[1][-1] != '*':
                                        if d[1] not in group_1 and d[1] not in group_2 and d[1] not in group_3 and d[1] not in group_4:
                                            group_4.append(d[1])
                                    elif right_star == False and left_star == False:
                                        if d[1] not in group_1 and d[1] not in group_2 and d[1] not in group_3 and d[1] not in group_4:
                                            group_4.append(d[1])
                                if d[1] == core_test:
                                    if right_star == True and d[0][0] != '*':
                                        if d[0] not in group_1 and d[0] not in group_2 and d[0] not in group_3 and d[0] not in group_4:
                                            group_4.append(d[0])
                                    elif left_star == True and d[0][-1] != '*':
                                        if d[0] not in group_1 and d[0] not in group_2 and d[0] not in group_3 and d[0] not in group_4:
                                            group_4.append(d[0])
                                    elif right_star == False and left_star == False:
                                        if d[0] not in group_1 and d[0] not in group_2 and d[0] not in group_3 and d[0] not in group_4:
                                            group_4.append(d[0])
                
                #if test_seq in whole_test:
                    #print('Group 4:', group_4)

                if len(group_1) > 0:
                    longest_core = get_longest(group_1, whole_seq)
                    core_whole_pairs.append([longest_core, whole_seq, fasta_seq, prots])
                    #print('Group 1 longest core:', longest_core)

                if len(group_2) > 0:
                    longest_core = get_longest(group_2, whole_seq)
                    core_whole_pairs.append([longest_core, whole_seq, fasta_seq, prots])
                    #print('Group 2 longest core:', longest_core)

                if len(group_3) > 0:
                    longest_core = get_longest(group_3, whole_seq)
                    core_whole_pairs.append([longest_core, whole_seq, fasta_seq, prots])
                    #print('Group 3 longest core:', longest_core)

                if len(group_4) > 0:
                    longest_core = get_longest(group_4, whole_seq)
                    core_whole_pairs.append([longest_core, whole_seq, fasta_seq, prots])
                    #print('Group 4 longest core:', longest_core)



    core_pairs_2 = []
    count = 0
    core_whole_pairs.sort(key=lambda x: x[2])
    fasta_ind = 0
    unique_wholes = []
    for a in core_whole_pairs:
        count += 1
        core = a[0]
        whole = a[1]
        fasta_seq = a[2]
        prot_str = a[3]
        unique_wholes.append(prot_str + '--' + whole + '--' + fasta_seq)

    unique_wholes = list(set(unique_wholes))
    unique_wholes.sort()

    for a in unique_wholes:
        prots = a.split('--')[0]
        whole = a.split('--')[1]
        fasta_seq = a.split('--')[2]

        core_matches = []
        for b in core_whole_pairs:
            core = b[0]
            whole_test = b[1]
            if whole_test == whole:
                core_test = re.sub('\*','',core)
                if core_test in whole and len(core_test) > 0:
                    print(whole, core)
                    start_pos = len(whole.split(core)[0])
                    end_pos = start_pos + len(core)
                    core_matches.append([core, start_pos, end_pos])

        bin_str = ''
        for i in range(0,len(whole)):
            pos_num = '0'
            for c in core_matches:
                if c[1] <= i <= end_pos:
                    pos_num = '1'
            bin_str += pos_num

        longest_cores = []
        ind_pos = 0
        for k in bin_str.split('0'):
            if len(k) == 0:
                ind_pos += 1
            else:
                longest_cores.append(whole[ind_pos:ind_pos+len(k)-1])
                ind_pos += len(k)

        for core in longest_cores:
            core_pairs_2.append([core, whole, fasta_seq, prots])


    unique_exps = []
    for a in peps:
        unique_exps.append(a[p_i[1]])
    unique_exps = list(set(unique_exps))
    unique_exps.sort()

    headers_out = ['Core Epitope', 'Proteins', 'Core Epitope Length', 'Whole Epitope', 'FASTA Seqs']
    for a in unique_exps:
        headers_out.append('% Rel. Intens. in ' + a)
    for a in unique_exps:
        headers_out.append('Intensity Landscape in ' + a)
    for a in unique_exps:
        headers_out.append('Unique Peptides in ' + a)

    out = []

    count = 0
    core_pairs_2.sort(key=lambda x: x[2])
    fasta_ind = 0
    unique_wholes = []
    for a in core_pairs_2:
        count += 1
        entry = ['']*len(headers_out)
        core = a[0]
        whole = a[1]
        fasta_seq = a[2]
        prot_str = a[3]

        entry[0] = core
        entry[1] = prot_str
        entry[2] = str(len(core))
        entry[3] = whole
        entry[4] = fasta_seq

        # 0 = experiment
        # 1 = pep
        # 2 = intensity
        pep_info = []
        core_pep_info = []

        for i in range(fasta_ind,len(peps)-1):
            b = peps[i]
            next_b = peps[i+1]
            if i == len(peps)-1:
                b = next_b
                exp = b[p_i[1]]
                peps_to_test = b[p_i[2]].split(';')
                intens_split = b[p_i[3]].split(';')
                for c in range(0,len(peps_to_test)):
                    pep_test = peps_to_test[c]
                    if pep_test != '' and pep_test in whole:
                        pep_info.append([exp, pep_test, intens_split[c]])
                        if core in pep_test or pep_test in core:
                            core_pep_info.append([exp, pep_test, intens_split[c]])
            if b[3] == fasta_seq:
                exp = b[p_i[1]]
                peps_to_test = b[p_i[2]].split(';')
                intens_split = b[p_i[3]].split(';')
                for c in range(0,len(peps_to_test)):
                    pep_test = peps_to_test[c]
                    if pep_test != '' and pep_test in whole:
                        pep_info.append([exp, pep_test, intens_split[c]])
                        if core in pep_test or pep_test in core:
                            core_pep_info.append([exp, pep_test, intens_split[c]])
                if next_b[3] != fasta_seq:
                    fasta_ind = i-1
                    break

        zeros = 0
        for exp in unique_exps:
            intens_ind = ''
            landscape_ind = ''
            pep_ind = ''
            for i in range(0,len(headers_out)):
                head_test = headers_out[i]
                if '% Rel. Intens. in ' + exp == head_test:
                    intens_ind = i
                if 'Intensity Landscape in ' + exp == head_test:
                    landscape_ind = i
                if 'Unique Peptides in ' + exp == head_test:
                    pep_ind = i

            pep_pairs = []
            for e in pep_info:
                if e[0] == exp:
                    pep_pairs.append([e[1], e[2]])

            core_pep_pairs = []
            for e in core_pep_info:
                if e[0] == exp:
                    core_pep_pairs.append([e[1], e[2]])
                    #print(core_pep_pairs[-1])

            core_intens = 0.0
            for e in core_pep_pairs:
                if e[1] != '':
                    core_intens += float(e[1])
            if core_intens == 0.0:
                zeros += 1

            whole_pos = []
            for i in range(0,len(whole)):
                whole_pos.append([i, whole[i]])

            pep_poses = []
            unique_peps = []
            for j in pep_pairs:
                pep_pos = []
                pep = j[0]
                if pep not in unique_peps:
                    unique_peps.append(pep)

                start_pos = len(whole.split(pep)[0])
                for i in range(start_pos,start_pos+len(pep)):
                    if j[1] != '':
                        pep_pos.append([i, float(j[1])])
                    else:
                        pep_pos.append([i, 0.0])
                pep_poses.append(pep_pos)

            intensity_landscape = []
            for i in range(0,len(whole_pos)):
                pos_intens = 0.0
                for j in pep_poses:
                    for k in j:
                        if k[0] == i:
                            pos_intens += float(k[1])
                intensity_landscape.append(pos_intens)
            landscape_str = ''
            for i in intensity_landscape:
                landscape_str += str(i) + '--'
            landscape_str = landscape_str[:-2]

            pep_str = ''
            for i in unique_peps:
                pep_str += i + ';'
            pep_str = pep_str[:-1]

            entry[intens_ind] = str(core_intens)
            entry[landscape_ind] = landscape_str
            entry[pep_ind] = pep_str

        if zeros != len(unique_exps):
            unique_wholes.append(whole)
            out.append(entry)
            print(entry[0], entry[1], count, '/', len(core_whole_pairs))


    unique_wholes = list(set(unique_wholes))
    unique_wholes.sort()
    out.sort(key=lambda x: x[3])
    out.append(['']*len(headers_out))

    out2 = out
    out = []
    start_ind = 0
    count = 0
    for whole in unique_wholes:
        count += 1
        cores = []
        prots = []
        fastas = []

        for i in range(start_ind,len(out2)-1):
            b = out2[i]
            next_b = out2[i+1]

            if b[3] == whole:
                cores.append(b[0])
                if b[1] not in prots:
                    prots.append(b[1])
                    fastas.append(b[4])
            
                if next_b[3] != whole:
                    start_ind = i-1
                    break
        cores2 = []
        for core in cores:
            if '*' not in core:
                if core+'*' in cores or '*'+core in cores:
                    pass
                else:
                    cores2.append(core)
            else:
                cores2.append(core)
        cores2 = list(set(cores2))

        prot_str = ''
        for prot in prots:
            prot_str += prot + ';'
        prot_str = prot_str[:-1]

        fasta_str = ''
        for f in fastas:
            fasta_str += f + ';'
        fasta_str = fasta_str[:-1]

        start_ind_2 = 0
        for core in cores2:
            entry = ['']*len(headers_out)
            entry[0] = core
            entry[1] = prot_str
            entry[2] = str(len(re.sub('\*','',core)))
            entry[3] = whole
            entry[4] = fasta_str

            for i in range(start_ind_2,len(out2)-1):
                b = out2[i]
                next_b = out2[i+1]

                if b[3] == whole:
                    if re.sub('\*','',b[0]) == re.sub('\*','',core):
                        for i in range(4,len(b)):
                            entry[i] = b[i]

                    if next_b[3] != whole:
                        start_ind_2 = i-1
                        break
            out.append(entry)
            print(entry[0], count, '/', len(unique_wholes))

    file_out = exp_name + '_core_epitopes_final.txt'
    savefile(file_out, out, headers_out)






def renorm(epi_file_in, imputation, filt_check):
    raw, headers, inds = openfile(epi_file_in, [])
    for i in range(0,len(headers)):
        if headers[i] == 'FASTA Seqs':
            last_ind = i + 1

    rel_cols = []
    for i in range(0,len(headers)):
        if '% Rel. Intens. in' in headers[i]:
            rel_cols.append(i)

    if filt_check == 'filt_ready' and imputation == 'no_imputation':
        raw2 = []
        for a in raw:
            zero_count = 0
            for i in range(last_ind,len(a)):
                if float(a[i]) == 0.0:
                    zero_count += 1
            if zero_count == 0:
                raw2.append(a)
        raw = raw2

    if imputation == 'lowest_all':
        all_vals = []
        for col in rel_cols:
            for a in raw:
                print(a)
                if float(a[col]) > 0.0:
                    all_vals.append(float(a[col]))
        if len(all_vals) > 0:
            impute_val = min(all_vals)
        else:
            impute_val = 1e-12

    out = []
    for a in raw:
        entry = ['']*len(headers)
        for i in range(0,len(a)):
            entry[i] = a[i]
        out.append(entry)

    count = 0


    #for col in range(last_ind,len(headers)):
    for col in rel_cols:
        rel_title = headers[col]
        exp = rel_title.split(' in ')[-1]
        landscape_col = ''
        for i in range(0,len(headers)):
            if 'Intensity Landscape in ' + exp in headers[i]:
                landscape_col = i

        count += 1
        all_vals = []
        for a in raw:
            if float(a[col]) > 0.0:
                all_vals.append(float(a[col]))

        total_sum = 0.0
        for a in raw:
            if float(a[col]) > 0.0:
                total_sum += float(a[col])
            else:
                if imputation == 'lowest_imputation':
                    total_sum += min(all_vals)
                elif imputation == 'lowest_all':
                    total_sum += impute_val

        new_sum = 0.0
        for i in range(0,len(raw)):
            val = raw[i][col]
            new_val = float(val)

            if new_val == 0.0 and imputation != 'no_imputation':
                if imputation == 'lowest_imputation':
                    new_val = min(all_vals)
                elif imputation == 'lowest_all':
                    new_val = impute_val

            if total_sum > 0.0:
                new_val = new_val / total_sum * 100.0
            
            new_sum += new_val
            out[i][col] = str(new_val)

            landscape = raw[i][landscape_col]
            new_landscape = ''
            for x in landscape.split('--'):
                new_val = float(x)
                if new_val == 0.0 and imputation != 'no_imputation':
                    if imputation == 'lowest_imputation':
                        new_val = min(all_vals)
                    elif imputation == 'lowest_all':
                        new_val = impute_val
                if total_sum > 0.0:
                    new_val = new_val / total_sum * 100.0
                new_landscape += str(new_val) + '--'
            new_landscape = new_landscape[:-2]
            out[i][landscape_col] = new_landscape
            
            print(new_val, count, '/', len(headers)+last_ind)

    out_file = epi_file_in.split('.txt')[0] + '_renorm.txt'

    if headers != []:
        out.insert(0,headers)

    filename = out_file.split('/')[-1]
    with open(os.path.expanduser(resultdir+filename), 'w') as a_file:
        for result in out:
            result = '\t'.join(result)
            a_file.write(result + '\n')



def final_output(exp_name, renorm_file, evidence_file, fasta_file, filt_check, dist_fig, venn_html, filt_params):
    cols = ['Core Epitope', 'Proteins', 'Core Epitope Length']
    filename = os.path.expanduser(resultdir+renorm_file.split('/')[-1])
    raw = list(csv.reader(open(filename, 'r'), delimiter='\t'))
    headers = raw[0]
    raw = raw[1:]

    count = 0
    inds = []
    for a in headers:
        for b in cols:
            if a == b:
                inds.append(count)

        count += 1

    cols = ['Sequence', 'Proteins']
    ev, e_h, e_i = openfile(evidence_file, cols)

    no_fasta_file = evidence_file.split('.txt')[0] + '_intens_norm_missing_fasta.txt'
    no_fasta, n_h, n_i = openfile(no_fasta_file, ['Missing'])

    unique_peps = []
    unique_prots = []
    for a in ev:
        unique_peps.append(a[e_i[0]])
        unique_prots.append(a[e_i[1]])
    unique_peps = list(set(unique_peps))
    unique_prots = list(set(unique_prots))


    to_print = []
    to_print.append('Contains <strong>' + str(len(unique_peps)) + '</strong> unique peptides from <strong>' + str(len(unique_prots)) + '</strong> unique proteins')

    fasta_file = str(fasta_file)
    fasta_file = re.sub("b'", '', fasta_file)
    fasta_file = re.sub("'", '', fasta_file)
    fasta_file = re.sub('_small', '', fasta_file)
    to_print.append('FASTA file used: ' + fasta_file.split('/')[-1])

    if filt_check == 'no_filt' or filt_check == 'test' or filt_check == 'quant_separately':
        to_print.append('No replicate filtering used.<br>')
    else:
        to_print.append('Filtering used:')
        to_print.append('Peptides must be quantified in at least '+str(filt_params[0])+' biological replicates.')
        to_print.append('Peptides must be quantified in at least '+str(filt_params[1])+' technical replicates.')
        if filt_params[2] == 'no_imputation':
            to_print.append('Missing values were not imputed.<br>')
        elif filt_params[2] == 'lowest_imputation':
            to_print.append('Missing values were imputed with the lowest measured value in each raw file.<br>')
        elif filt_params[2] == 'lowest_all':
            to_print.append('Missing values were imputed with the lowest measured value in all combined raw files.<br>')

    no_fasta_str = ''
    for x in no_fasta:
        if len(x) > 0:
            no_fasta_str += x[0] + ', '
    no_fasta_str = no_fasta_str[:-2]

    if len(no_fasta) > 0:
        to_print.append('<small>NB: peptides from the following <b>'+str(len(no_fasta))+'</b> proteins were not found in the given FASTA file, and have therefore been omitted from the analysis:<br> <i>' + no_fasta_str + '</i></small><br>')

        to_print.append('<br>')

    unique_prots = []
    lens = []
    for a in raw:
        unique_prots.append(a[inds[1]])
        lens.append(float(a[inds[2]]))
    unique_prots = list(set(unique_prots))
    len_sum = 0.0
    for b in lens:
        len_sum += b
    if len(lens) > 0:
        avg_len = len_sum / float(len(lens))
    else:
        avg_len = 0.0

    to_print.append('<b>PLAtEAU Analysis</b><br><strong>'+str(len(raw)) + '</strong> unique core epitopes found from <strong>' + str(len(unique_prots)) + '</strong> unique proteins')
    to_print.append('Average core epitope length: ' + str(avg_len)[:6] + ' residues')

    to_print.append('<br>')
    if dist_fig != 'no_peptides':
        dist_fig = resultdir + dist_fig.split('/')[-1]
        to_print.append('Distribution of unique peptide and epitope lengths:')
        to_print.append('<img src="'+dist_fig+'"><br>')
    else:
        to_print.append('No peptides found matching the specified criteria!<br>')

    to_print.append('<br>')

    to_print.append('Download output file <a href="'+resultdir+renorm_file.split('/')[-1]+'" download="'+renorm_file+'" target="_blank" rel="noopener">here</a>.')

    if len(venn_html) > 0:
        for a in venn_html:
            to_print.append(a)
 
    html_output(to_print)	



def get_exps(evidence_file, params, template_file):
    cols = ['Raw file']
    raw, headers, inds = openfile(evidence_file, cols)

    if len(template_file.split('/')[-1]) > 0:
        template, t_h, t_i = openfile(template_file, ['Raw File'])
    else:
        template = ['lll','lll','lll','lll']

    unique_exps = []
    for a in raw:
        if a[inds[0]] not in unique_exps:
            unique_exps.append(a[inds[0]])

    unique_exps.sort()
    exp_matches = []
    for a in unique_exps:
        exp_cond = 'Condition'
        exp_bio = 'Bio. rep.'
        exp_tech = 'Tech. rep.'

        for b in template:
            if len(b) > 0:
                if b[0] == a:
                    exp_cond = b[1]
                    exp_bio = b[2]
                    exp_tech = b[3]
        exp_matches.append([a, exp_cond, exp_bio, exp_tech])
	
    to_print = []
    to_print.append('<center>')
    to_print.append('<div class="text-box-upload-and-run">')
    to_print.append('<form action="../php/run-script_2.php" method="post" enctype="multipart/form-data">')
    to_print.append('<input type="hidden" name="expdir" value="'+expdir+'"/>')

    exp_lens = []
    for a in unique_exps:
        exp_lens.append(len(a))
    if len(exp_lens) > 0:
        place_name = '&nbsp;'*max(exp_lens)
    else:
        place_name = '&nbsp;'
    spacer = '&nbsp;'*8

    vars = []

    if len(unique_exps) > 0:
        to_print.append('Fill out the organization of biological and technical replicates.<br><br><i>Example:</i><br> <img src="../images/example_table.png" width=400>')

        to_print.append('<font face="Courier">')

    for i in range(1,len(unique_exps)+1):
        exp = unique_exps[i-1]
        cond_name = 'cond_' + str(i)
        bio_rep_name = 'bio_rep_' + str(i)
        tech_rep_name = 'tech_rep_' + str(i)

        cond_val = exp_matches[i-1][1]
        bio_val = exp_matches[i-1][2]
        tech_val = exp_matches[i-1][3]

        cond_form = '<input type="text" name="'+cond_name+'" id="'+cond_name+'" style="width:20%" value="'+cond_val+'" onfocus="if (this.value==\''+cond_val+'\') this.value=\'\';"/>'
        bio_form = '<input type="text" name="'+bio_rep_name+'" id="'+bio_rep_name+'" style="width:10%" value="'+bio_val+'" onfocus="if (this.value==\''+bio_val+'\') this.value=\'\';"/>'
        tech_form = '<input type="text" name="'+tech_rep_name+'" id="'+tech_rep_name+'" style="width:10%" value="'+tech_val+'" onfocus="if (this.value==\''+tech_val+'\') this.value=\'\';"/>'

        to_print.append(exp + '&emsp; ' + cond_form + ' ' + bio_form + ' ' + tech_form)		

        vars.append([exp+';'+cond_name+';'+bio_rep_name+';'+tech_rep_name+';x'])

    to_print.append('</font>')
    to_print.append('<br>')

    param_file = params[0][0].split(';')[0] + '_params.txt'
    savefile(os.path.expanduser(expdir+param_file), params, [])
    var_file = params[0][0].split(';')[0] + '_vars.txt'
    savefile(os.path.expanduser(expdir+var_file), vars, [])

    # Imputation
    if len(unique_exps) > 0:
        to_print.append('Imputation of missing values. <i>NB: no imputation means only peptides found in <b>all</b> biological and technical replicates will be included</i><br>')

        to_print.append('<input type="radio" name="imputation" value="no_imputation">No imputation.')
        to_print.append('<input type="radio" name="imputation" value="lowest_imputation" checked>Impute with lowest measured value in each run.')
        to_print.append('<input type="radio" name="imputation" value="lowest_all">Impute with lowest measured value of all combined runs.<br>')

	
        to_print.append('<input type="hidden" name="param_file" value="'+param_file+'"/>')
        to_print.append('<input type="hidden" name="var_file" value="'+var_file+'"/>')

        to_print.append('<input type="submit" value="Next Page" name="submit">')
    else:
        to_print.append('It seems like the MaxQuant evidence file is not formatted properly, or is empty. Please check it.<br>')


    html_output(to_print)



def gen_len_dist(exp, evidence_file, epitope_file):
    cols = ['Sequence']
    evidence, e_h, e_i = openfile(os.path.expanduser(evidence_file), cols)

    cols = ['Core Epitope']
    epitope_file = os.path.expanduser(resultdir+epitope_file.split('/')[-1])
    epitopes, ep_h, ep_i = openfile(os.path.expanduser(epitope_file), cols)

    unique_peptides = []
    peptide_lens = []
    for a in evidence:
        unique_peptides.append(a[e_i[0]])	
    unique_peptides = list(set(unique_peptides))
    for a in unique_peptides:
        peptide_lens.append(len(a))

    unique_epitopes = []
    epitope_lens = []
    for a in epitopes:
        unique_epitopes.append(a[ep_i[0]])	
    unique_epitopes = list(set(unique_epitopes))
    for a in unique_epitopes:
        if '*' in a:
            a = re.sub('\*','',a)
        epitope_lens.append(len(a))

    fig = plt.figure()
    ax = plt.subplot(111)

    if len(peptide_lens) > 0:
        max_len = max(peptide_lens)
        min_len = min(peptide_lens)
        num_bins = max_len - min_len
    
    print('Bins:', num_bins, len(peptide_lens), len(epitope_lens))

    if len(peptide_lens) > 1 and num_bins > 0:
        n, bins, patches = plt.hist(peptide_lens, max(peptide_lens)-min(peptide_lens), facecolor='blue', alpha=0.5, label='Peptides')
        if len(epitope_lens) > 1 and num_bins > 0:
            n, bins, patches = plt.hist(epitope_lens, max(epitope_lens)-min(epitope_lens), facecolor='red', alpha=0.5, label='Epitopes')
	
        ax.set_xlabel('Length (AA)')
        ax.set_ylabel('Number')
        plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
	
        out_file = os.path.expanduser(resultdir+exp+'_pep_epi_dist.svg')
        plt.savefig(out_file, bbox_inches='tight', transparent=True)
    else:
        out_file = 'no_peptides'

    return out_file




def venn_2(set_name, name_a, name_b, list_a, list_b, exp):
    a_num = str(len(list_a))
    b_num = str(len(list_b))
		
    a_set = set(list_a)
    b_set = set(list_b)

    fig_x = 2.5
    fig_y = 2.5

    plt.rcParams.update({'font.size': 8})
    plt.figure(figsize=(fig_x,fig_y))
    v = venn2([a_set, b_set], (name_a + ' (' + a_num + ')', name_b + ' (' + b_num + ')'))
    for text in v.set_labels:
        text.set_fontsize(8)	
	
    plt.rcParams.update({'font.size': 12})
    plt.title(set_name)

    out_file = resultdir + exp + '_' + set_name + '_' + name_a + '_' + name_b + '_venn.svg'
    plt.savefig(os.path.expanduser(out_file), transparent=True)
    plt.clf()

    return out_file


def venn_2_conds(exp, cond1, cond2, var_raw, bio_rep_min, tech_rep_min, pep_pass_file, evidence_file, renorm_file):
    vars = []
    for a in var_raw:
        raw_file = a[0]
        cond = a[1]
        bio_rep = a[2]
        tech_rep = a[3]

    pep_pass_table, pep_headers, p_i = openfile(os.path.expanduser(pep_pass_file), [])

    html_out = []
    html_out.append('<br><hr><b>' + cond1 + ' vs. ' + cond2 + '</b>')

    cols = ['Sequence', 'Proteins']
    evidence, e_h, e_i = openfile(os.path.expanduser(evidence_file), cols)

    cols = ['Core Epitope', 'Proteins', 'Whole Epitope']
    renorm, r_h, r_i = openfile(os.path.expanduser(renorm_file), cols)

    # Unique passing peptides
    for i in range(0,len(pep_headers)):
        if 'Pass in ' + cond1 in pep_headers[i]:
            ind1 = i
        if 'Pass in ' + cond2 in pep_headers[i]:
            ind2 = i
    peps_1 = []
    peps_2 = []
    for a in pep_pass_table:
        if a[ind1] == 'Yes':
            peps_1.append(a[0])
        if a[ind2] == 'Yes':
            peps_2.append(a[0])

    pep_fig_file = venn_2('Unique passing peptides', cond1, cond2, peps_1, peps_2, exp)
    html_out.append('<img src="'+pep_fig_file+'"')

    # Proteins
    prots_1 = []
    prots_2 = []
    for a in evidence:
        if a[e_i[0]] in peps_1:
            if a[e_i[1]] not in prots_1:
                prots_1.append(a[e_i[1]])
        if a[e_i[0]] in peps_2:
            if a[e_i[1]] not in prots_2:
                prots_2.append(a[e_i[1]])

    prot_fig_file = venn_2('Unique passing proteins', cond1, cond2, prots_1, prots_2, exp)
    html_out.append('<img src="'+prot_fig_file+'"')

    # Core epitopes
    epis_1 = []
    epis_2 = []
    for a in renorm:
        whole = a[r_i[2]]
        for b in peps_1:
            if b in whole and a[r_i[0]] not in epis_1:
                epis_1.append(a[r_i[0]])
        for b in peps_2:
            if b in whole and a[r_i[0]] not in epis_2:
                epis_2.append(a[r_i[0]])

    epi_fig_file = venn_2('Core epitopes', cond1, cond2, epis_1, epis_2, exp)
    html_out.append('<img src="'+epi_fig_file+'"')
	
    return html_out



def volcano_2_conds(exp, cond1, cond2, var_raw, bio_rep_min, tech_rep_min, pep_pass_file, evidence_file, renorm_file, ratio_cutoff, p_cutoff):
    html_out = []

    vars = []
    raws1 = []
    raws2 = []

    for a in var_raw:
        raw_file = a[0]
        cond = a[1]
        bio_rep = a[2]
        tech_rep = a[3]
        if cond == cond1:
            raws1.append(raw_file)
        if cond == cond2:
            raws2.append(raw_file)

    pep_pass_table, pep_headers, p_i = openfile(os.path.expanduser(pep_pass_file), [])

    cols = ['Core Epitope', 'Proteins', 'Whole Epitope']
    renorm, r_h, r_i = openfile(os.path.expanduser(renorm_file), cols)

    inds1 = []
    inds2 = []
    for i in range(0,len(r_h)):
        for a in raws1:
            if '% Rel. Intens. in ' + a in r_h[i]:
                inds1.append(i)
        for a in raws2:
            if '% Rel. Intens. in ' + a in r_h[i]:
                inds2.append(i)

    out = []
    headers_out = []
    for i in range(0,4):
        headers_out.append(r_h[i])

    for i in range(4,len(r_h)):
        if i in inds1 or i in inds2:
            #headers_out.append('log2('+r_h[i]+')')
            headers_out.append(r_h[i])
	
    headers_out.append('log2((' + cond1 + ')/(' + cond2 + ') mean ratio)')
    headers_out.append("'-log10((" + cond1 + ')/(' + cond2 + ') p)')

    renorm2 = []
    for a in renorm:
        entry = ['']*len(headers_out)
        for i in range(0,4):
            entry[i] = a[i]

        for i in range(0,len(inds1)):
            entry[4+i] = a[inds1[i]]
        for i in range(0,len(inds2)):
            entry[4+len(inds1)+i] = a[inds2[i]]
        renorm2.append(entry)

    inds1 = []
    inds2 = []
    for i in range(0,len(headers_out)):
        for a in raws1:
            if '% Rel. Intens. in ' + a in headers_out[i]:
                inds1.append(i)
        for a in raws2:
            if '% Rel. Intens. in ' + a in headers_out[i]:
                inds2.append(i)
    ratio_means = []
    p_vals = []
    upreg_pass = []
    downreg_pass = []
    upreg_means = []
    upreg_ps = []
    downreg_means = []
    downreg_ps = []
    abs_ratios = []

    for a in renorm2:
        entry = ['']*len(headers_out)
        for i in range(0,len(a)):
            entry[i] = a[i]

        vals1 = []
        vals2 = []
        ratios = []
        ratios_log2 = []

        for i in inds1:
            #val = math.log(float(a[i]),2)
            val = float(a[i])
            vals1.append(val)
            #entry[i] = str(math.log(val,2))
            entry[i] = str(val)

        for i in inds2:
            #val = math.log(float(a[i]),2)
            val = float(a[i])
            vals2.append(val)
            #entry[i] = str(math.log(val,2))
            entry[i] = str(val)

        for i in range(0,len(vals1)):
            ratio = vals1[i] / vals2[i]
            ratios.append(ratio)
	
            ratio = math.log(vals1[i]/vals2[i],2)
            ratios_log2.append(ratio)
            abs_ratios.append(abs(ratio))

        if len(ratios) > 0:
            mean = sum(ratios)/float(len(ratios))
        else:
            mean = 0.0
        if len(ratios_log2) > 0:
            mean_log2 = sum(ratios_log2)/float(len(ratios_log2))
        else:
            mean_log2 = 0.0
        p = stats.ttest_1samp(ratios_log2, 0)[1]
        p_log10 = -math.log(float(p),10)

        if p <= p_cutoff:
            if mean >= ratio_cutoff:
                upreg_pass.append([mean_log2, p_log10])
                upreg_means.append(mean_log2)
                upreg_ps.append(p_log10)

            if mean <= 1.0/ratio_cutoff:
                downreg_pass.append([mean_log2, p_log10])
                downreg_means.append(mean_log2)
                downreg_ps.append(p_log10)

        ratio_means.append(mean_log2)
        p_vals.append(p_log10)
        entry[-2] = str(mean_log2)
        entry[-1] = str(p)
        out.append(entry)

    table_out_file = resultdir + exp + '_' + cond1 + '_' + cond2 + '_volcano_table.txt'
    savefile(os.path.expanduser(table_out_file), out, headers_out)


    # Volcano plot
    plt.rcParams.update({'font.size': 20})
	
    fig_x = 10
    fig_y = 10
    margin = 5

    #y_min = min(p_vals) / 10.0
    #y_max = max(p_vals) * 10.0

    if len(abs_ratios) > 0:
        x_max = max(abs_ratios) + 2.0
    else:
        x_max = 10.0
    x_min = -x_max

    x_cutoff = -math.log(p_cutoff, 10)
    y_cutoff_1 = math.log(ratio_cutoff,2)
    y_cutoff_2 = math.log((1.0/ratio_cutoff),2)

    fig = plt.figure(figsize=(fig_x,fig_y))
    ax = plt.gca()

    plt.axvline(x=y_cutoff_1, color="red", linewidth=3)
    plt.axvline(x=y_cutoff_2, color="red", linewidth=3)
    plt.axhline(y=x_cutoff, color="red", linewidth=3)

    #rect_start = (y_cutoff, x_cutoff)
    #rect_width = x_max - y_cutoff
    #rect_height = y_max - x_cutoff
    #ax.add_patch(patches.Rectangle(rect_start, rect_width, rect_height, color="red", alpha=0.1, zorder=-1))

    ax.scatter(ratio_means, p_vals, color="grey", s=30)
    ax.scatter(upreg_means, upreg_ps, color="blue", s=30)
    ax.scatter(downreg_means, downreg_ps, color="blue", s=30)

    #ax.set_yscale('log')

    #ax.set_ylim([y_min, y_max])
    ax.set_xlim([x_min, x_max])
	
    #plt.axvline(x=x_min, color='black', linewidth=4)
    #plt.axvline(x=x_max, color='black', linewidth=4)
    #plt.axhline(y=y_min, color='black', linewidth=4)
    #plt.axhline(y=y_max, color='black', linewidth=4)

    plt.xlabel('log2((' + cond1 + ')/(' + cond2 + ') mean ratio)')
    plt.ylabel('-log10(p value)')

    #plt.text(1, 1, str(len(upreg_means)), color="blue", fontsize=20)
    #plt.text(0, 1, str(len(downreg_means)), color="blue", fontsize=20)
    ax.annotate(str(len(upreg_means)), xy=(0.90,0.95), xycoords='axes fraction', color="blue", fontsize=20)
    ax.annotate(str(len(downreg_means)), xy=(0.05,0.95), xycoords='axes fraction', color="blue", fontsize=20)

    fig_file_out = resultdir + exp + '_' + cond1 + '_' + cond2 + '_volcano_plot.svg'
    plt.savefig(os.path.expanduser(fig_file_out), transparent=True)

    html_out.append('<img src="'+fig_file_out+'" width=600>')
    html_out.append('<a href="'+resultdir+table_out_file.split('/')[-1]+'" download="'+table_out_file.split(resultdir)[-1]+'" target="_blank" rel="noopener">Download volcano plot table</a>')
	
    return html_out




def check_file(file_in, banned_chars, banned_files):
    check = True
    name_check = file_in.split('/')[-1]

    if len(name_check) == 0:
        check = False

    for x in name_check:
        if x in banned_chars:
            check = False

    for x in banned_files:
        if x in file_in:
            check = False

    return check





# Regenerate the whole landscape for each protein, to include jumps 
def remodel_wholes(final_file, fasta_file, epi_file):
    cols = ['Core Epitope', 'Proteins', 'Whole Epitope']
    raw, headers, inds = openfile(final_file, cols)

    cols = ['Protein', 'Experiment', 'FASTA seq', 'Core epitopes', 'Whole epitopes', 'PLAtEAU Intens. (norm)', 
            'Peptides contributing to core epitopes']
    cores, c_h, c_i = openfile(epi_file, cols)

    unique_proteins = []
    for a in raw:
        for b in a[inds[1]].split(';'):
            unique_proteins.append(b)
    unique_proteins = list(set(unique_proteins))

    prot_fasta = []
    unique_peptides = []
    for a in cores:
        for b in a[c_i[6]].split(';'):
            for c in b.split(','):
                if c != '':
                    unique_peptides.append(c)
        prot_fasta.append(a[c_i[0]] + '--' + a[c_i[2]])
    unique_peptides = list(set(unique_peptides))
    prot_fasta = list(set(prot_fasta))
    
    unique_proteins.sort()
    prot_fasta.sort()
    start_ind = 0

    new_wholes = []

    count = 0
    for prot in unique_proteins:
        count += 1
        fasta_seq = ''
        for i in range(start_ind,len(prot_fasta)):
            split1 = prot_fasta[i].split('--')
            if split1[0] == prot:
                fasta_seq = split1[1]
                start_ind = i-1
                break

        prot_seq = []
        # 0 = pos, 1 = res, 2 = peps found
        for i in range(0,len(fasta_seq)):
            prot_seq.append([i, fasta_seq[i], 0])
        pos_to_add = []

        for pep in unique_peptides:
            if pep in fasta_seq:
                start_pos = len(fasta_seq.split(pep)[0])
                end_pos = start_pos + len(pep)
                #print(pep)
                #print(fasta_seq[start_pos:end_pos])
                #print('')
                for i in range(start_pos,end_pos):
                    pos_to_add.append(i)

        prot_seq_2 = []
        reses = []
        nums = []
        res_str = ''
        for i in range(0,len(prot_seq)):
            pos = prot_seq[i][0]
            res = prot_seq[i][1]
            num = 0
            for j in pos_to_add:
                if pos == j:
                    num += 1
            prot_seq_2.append([pos, res, num])
            reses.append(res)
            nums.append(num)
            res_str += res
            #print(prot_seq_2[-1])

        res_str_2 = ''
        for i in range(0,len(reses)):
            if nums[i] != 0:
                res_str_2 += reses[i]
            else:
                res_str_2 += '-'

        whole_split = res_str_2.split('-')
        wholes = []
        for a in whole_split:
            if a != '':
                wholes.append(a)
                #print(a, count, '/', len(unique_proteins))
                new_wholes.append(a)

    new_wholes = list(set(new_wholes))
    final_out = []

    for a in raw:
        entry = ['']*len(a)
        for i in range(0,len(a)):
            entry[i] = a[i]

        core = a[inds[0]]
        core_test = re.sub('\*', '', core)
        #print('')
        #print('Core:', core)

        for b in new_wholes:
            if core_test in b:
                entry[inds[2]] = b
                #print('Whole:', b)

        final_out.append(entry)

    out_file = final_file.split('.txt')[0] + '_refit.txt'
    savefile(out_file, final_out, headers)








def iedb_lookup(pep, dr_allele, iedb_path):
    if dr_allele != 'DR':
        first = dr_allele.split('_')[0]
        second = dr_allele.split('_')[-1]
        alpha = second[:2]
        beta = second[2:]
        test_1 = first + '*' + alpha + ':' + beta
    else:
        test_1 = 'DR'

    if iedb_path[-1] != '/':
        iedb_path = iedb_path + '/'

    files_to_open = []
    file_list = os.listdir(iedb_path)
    for a in file_list:
        if test_1 in a:
            files_to_open.append(os.path.expanduser(iedb_path + a))
    
    raw_comb = []
    for iedb_file in files_to_open:
        raw = list(csv.reader(open(iedb_file, 'r'), delimiter='\t'))
        if len(raw) > 0:
            headers = raw[0]
            raw = raw[1:]
        else:
            headers = []
            raw = []

        new_headers = []
        for i in range(0,len(headers)):
            if headers[i] not in new_headers:
                new_headers.append(headers[i])
            else:
                new_headers.append(headers[i] + '-2')

        cols = [
            'PubMed ID', # 0
            'Authors', # 1
            'Journal', # 2
            'Date', # 3
            'Title', # 4
            'Epitope IRI', # 5 IEDB link
            'Description', # 6 Peptide sequence
            'Units', # 7
            'Qualitative Measure', # 8
            'Measurement Inequality', # 9  > etc
            'Quantitative measurement', # 10
            'Allele Name', # 11
            'MHC allele class', # 12
            ]


        count = 0
        inds = []
        for a in new_headers:
            for b in cols:
                if a == b:
                    inds.append(count)
        count += 1

        for a in raw:
            raw_comb.append(a)

    out = []

    for a in raw_comb:
        if len(a) > inds[-1]:
            if pep == a[inds[6]]:
                units = a[inds[7]]
                qual = a[inds[8]]
                inequal = a[inds[9]]
                quant = a[inds[10]]
                if quant != '':
                    quant = quant + ' ' + units
                if inequal != '':
                    quant = inequal + ' ' + quant

                authors = a[inds[1]]
                pmid = a[inds[0]]
                year = a[inds[3]]

                authors = authors.split(';')[0]
                authors = authors.split(' ')[-1]
                cite = authors + ' et al. (' + year + ') (PMID: ' + pmid + ')'
                    
                iedb_num = a[inds[5]]

                out = [pep, a[inds[9]], cite, iedb_num, qual, quant]

    return out




def draw_pep(pep_seq, pep_color, label_text, img_width, res_spacing):
    #fontname_seq = 'monaco.ttf' # for seq
    fontname_seq = 'DejaVuSans.ttf'
    fontsize_seq = 85
    fonttyp_seq = ImageFont.truetype(fontname_seq, fontsize_seq)

    fontname_label = 'Calibri.ttf' # for label
    fontsize_label = 40
    fonttyp_label = ImageFont.truetype(fontname_label, fontsize_label)

    pep_seq = re.sub(' ', '', pep_seq)

    background = (255,255,255)

    right_side = 800
    #right_side = 20
    scrap_width = (100*len(pep_seq)) + 8000

    image = Image.new("RGB", (scrap_width,500), background)
    draw = ImageDraw.Draw(image)
    margin_top = 3 * 2
    margin_side = 20 * 2
    w, h = draw.textsize(pep_seq, font=fonttyp_seq)
    w += (res_spacing*len(pep_seq))

#    if img_width == 18:
#        width = 1608
#    else:
#        width = w + margin_side
    
    width = 3000
    height = h + margin_top

    image2 = Image.new("RGB", (width+right_side, height+8000), background)

    fill = "   "
    x = 0
    w_fill, y = draw.textsize(fill, font=fonttyp_seq)
    y_mid = (height-h)/2
    #x_draw, x_paste = 20, 300
    x_draw, x_paste = 20, 20


    for i in pep_seq:
        w_full = draw.textsize(fill+i, font=fonttyp_seq)[0]
        w_full += res_spacing
        w = w_full - w_fill # width of character

        if i != ' ':
            draw.rectangle((x_draw+w_fill, 0, x_draw+w_full, y+5), pep_color)
            draw.text(((x_draw+int(float(res_spacing)/2.0)),-5), fill+i, (0,0,0), font=fonttyp_seq)
        else:
            draw.rectangle((x_draw+w_fill, 0, x_draw+w_full, y+5), background)
            draw.text((x_draw,0), fill+i, background, font=fonttyp_seq)

        iletter = image.crop((x_draw+w_fill, 0, x_draw+w_full, y+5))
        image2.paste(iletter, (x_paste, 0))
        x_draw += w_full
        x_paste += w

    draw2 = ImageDraw.Draw(image2)
    #draw2.text((width+300, 0), label_text, (0,0,0), font=fonttyp_label)


    cropped = image2.crop((0,0,width+right_side,height))
    out_name = os.path.expanduser('tmp/'+label_text + '_' + pep_seq + '_peptide.png')
    out_name = re.sub(' ', '', out_name)
    cropped.save(out_name)

    return out_name






def netmhciipan(final_file_in, length, dr_allele, core_file, fasta_file_whole, high_cutoff, med_cutoff, temp_path, net_path, epi_file, iedb_path):
    cols = ['Core Epitope', 'Proteins', 'Whole Epitope', 'FASTA Seqs']
    raw, headers, inds = openfile(os.path.expanduser(final_file_in), cols)

    cols = ['Protein', 'Experiment', 'FASTA seq', 'Core epitopes', 'Whole epitopes', 'PLAtEAU Intens. (norm)', 
            'Peptides contributing to core epitopes']
    cores, c_h, c_i = openfile(epi_file, cols)

    unique_exps = []
    for i in headers:
        if '% Rel. Intens. in ' in i:
            exp = i.split('% Rel. Intens. in ')[-1]
            unique_exps.append(exp)
    unique_exps = list(set(unique_exps))

    raw_split = []
    for a in raw:
        if ';' in a[inds[1]] and ';' in a[inds[3]]:
            prot_split = a[inds[1]].split(';')
            fasta_split = a[inds[3]].split(';')
            print(prot_split)
            print(fasta_split)
            for i in range(0,len(prot_split)):
                entry = ['']*len(a)
                for j in range(0,len(a)):
                    entry[j] = a[j]
                entry[inds[1]] = prot_split[i]
                entry[inds[3]] = fasta_split[i]
                raw_split.append(entry)
        else:
            raw_split.append(a)

    raw_split.sort(key=lambda x: x[4])
    
    prot_fasta = []
    for a in cores:
        prot_fasta.append(a[c_i[0]] + '--' + a[c_i[2]])
    prot_fasta = list(set(prot_fasta))
    
    lens = []
    for a in raw:
        lens.append(len(a[inds[0]]))
    if len(lens) > 0:
        longest_len = max(lens)
    else:
        longest_len = 16

    unique_wholes = []

    extra_side = 8 # how many extra residues on each side

    # all sequences should be this length
    master_len = extra_side + longest_len + extra_side
    min_len = 28
    if master_len < min_len:
        master_len = min_len
        longest_len = longest_len - (extra_side*2)

    count = 0
    for a in raw_split:
        count += 1
        entry = ['']*len(headers)
        for i in range(0,len(a)):
            entry[i] = a[i]

        whole_seq = a[inds[2]]
        whole_seq = re.sub('\*', '', whole_seq)

        fig_margin_right = extra_side
        fig_margin_left = extra_side

        margin_spacer_right = '_' * fig_margin_right
        margin_spacer_left = '_' * fig_margin_left
        
        unique_whole_seqs = []
        prot_fasta_rev = []

        sequence = a[inds[3]]
        seq_start = sequence.split(whole_seq)[0][-fig_margin_left:]
        seq_end = sequence.split(whole_seq)[-1][:fig_margin_right]
        whole_seq_with_margin = seq_start + '-' + whole_seq + '-' + seq_end
        unique_whole_seqs.append(whole_seq_with_margin)

        prot_str = a[inds[1]]      

        #for b in prot_fasta:
        #    split1 = b.split('--')
        #    sequence = split1[1]
        #    prot_fasta_rev.append(sequence + '--' + split1[0])
        #    if whole_seq in sequence:
        #        seq_start = sequence.split(whole_seq)[0][-fig_margin_left:]
        #        seq_end = sequence.split(whole_seq)[-1][:fig_margin_right]
        #        whole_seq_with_margin = seq_start + '-' + whole_seq + '-' + seq_end
        #        unique_whole_seqs.append(whole_seq_with_margin)
        #unique_whole_seqs = list(set(unique_whole_seqs))
        #unique_whole_seqs.sort()
        #prot_fasta_rev.sort()
        #print(a[0], len(unique_whole_seqs), count, '/', len(raw))

        #start_ind = 0
        for b in unique_whole_seqs:
            seq_start = b.split('-')[0]
            seq_end = b.split('-')[-1]
            seq_center = b.split('-')[1]
        #    prot_matches = []
        #    for i in range(start_ind,len(prot_fasta_rev)):
        #        split1 = prot_fasta_rev[i].split('--')
        #        if seq_start+seq_center+seq_end in split1[0]:
        #            prot_matches.append(split1[1])
        #            start_ind = i-1
        #            break
        #    prot_matches = list(set(prot_matches))
        #    prot_str = ''
        #    for x in prot_matches:
        #        prot_str += x + ';'
        #    prot_str = prot_str[:-1]

            if len(seq_start) < fig_margin_left:
                delt_val = fig_margin_left - len(seq_start)
                add_marg = margin_spacer_left[:delt_val]
                seq_start = add_marg + seq_start
            if len(seq_end) < fig_margin_right:
                delt_val = fig_margin_right - len(seq_end)
                add_marg = margin_spacer_right[:delt_val]
                seq_end = seq_end + add_marg

            whole_seq = seq_start + seq_center + seq_end
            unique_wholes.append(sequence + '--' + prot_str + '--' + whole_seq + '--' + seq_start + '--' + seq_end)
            #print(unique_wholes[-1])

    unique_wholes = list(set(unique_wholes))
    unique_wholes.sort()
    raw_split.append(['']*len(headers))

    count = 0
    fasta_ind = 0
    for x in unique_wholes:
        count += 1
        b = x.split('--')
        fasta_seq = b[0]
        prot = b[1]
        whole_seq = b[2]
        seq_start = b[3]
        seq_end = b[4]
        cores_b = []
        passing_exps = []

        for i in range(fasta_ind,len(raw_split)-1):
            b = raw_split[i]
            next_b = raw_split[i+1]

            if b[inds[3]] == fasta_seq:
                core_epi = b[0]
                if '*' in core_epi:
                    core_epi = re.sub('\*','',core_epi)
                if core_epi in whole_seq:
                    cores_b.append(core_epi)
                    for i in range(0,len(a)):
                        if 'Unique Peptides in ' in headers[i] and a[i] != '':
                            exp = headers[i].split('Unique Peptides in ')[-1]
                            passing_exps.append(exp)

                if next_b[inds[3]] != fasta_seq:
                    fasta_ind = i-1
                    break

        #for a in raw_split:
        #    core_epi = a[0]
        #    core_epi = re.sub('\*', '', core_epi)
        #    if core_epi in whole_seq and core_epi not in cores_b:
        #        cores_b.append(core_epi)
        #        for i in range(0,len(a)):
        #            if 'Unique Peptides in ' in headers[i] and a[i] != '':
        #                exp = headers[i].split('Unique Peptides in ')[-1]
        #                passing_exps.append(exp)

        cores_b = list(set(cores_b))
        core_tables = []
        core_poses = []
        for core_epi in cores_b:
            first_pos = len(whole_seq.split(core_epi)[0])
            last_pos = first_pos + len(core_epi)
            core_table = []
            core_pos = []
            for i in range(first_pos,last_pos):
                core_table.append([whole_seq[i], i])
                core_pos.append(i)
            core_tables.append(core_table)
            core_poses.append(core_pos)


        for exp in passing_exps:
            net_test_seq = whole_seq
            if len(seq_start) > 0:
                net_test_seq = net_test_seq.split(seq_start)[-1]
            if len(seq_end) > 0:
                net_test_seq = net_test_seq.split(seq_end)[0]
            net_test_seq = net_test_seq[1:-1]
            fasta_out = [['>'+net_test_seq], [net_test_seq]]
            net_name = 'nettemp'
            fasta_file = temp_path + net_name + '.fsa'
            savefile(os.path.expanduser(fasta_file), fasta_out, [])

            net_file = net_name + '_save_in_specified.fs.out'

            net_prog = os.path.expanduser(net_path + 'netMHCIIpan')
            in_file = os.path.expanduser(fasta_file)
            out_file = os.path.expanduser(temp_path + net_file)

            cmd = net_prog + ' -f ' + in_file + ' -a ' + dr_allele + ' -length ' + str(length) + ' >| ' + out_file
            os.system(cmd)
            net_table, n_h, n_i = openfile(os.path.expanduser(temp_path + net_file), [])

            inf_table = [] # 0 = whole_seq, 1 = pep, 2 = 1-log50(aff), 3 = pos_str

            for c in net_table:
                if len(c) > 0:
                    split_str = c[0].split(' ')
                    new_a = []
                    if len(split_str) > 0:
                        for b in split_str:
                            if b != '':
                                new_a.append(b)

                        if new_a[0] == 'Seq':
                            #cols = ['Peptide', '1-log50k(aff)']
                            cols = ['Peptide', 'Affinity(nM)', '%Rank']
                            inds_2 = []
                            for j in cols:
                                for i in range(0,len(new_a)):
                                    if j == new_a[i]:
                                        inds_2.append(i)

                        if len(new_a) > 2:
                            if new_a[1] == dr_allele:
                                pep = new_a[inds_2[0]]
                                #rev_inf = new_a[inds_2[1]]
                                rev_inf = 1.0/(float(new_a[inds_2[1]])) # Ka, inverse of Kd

                                first_pos = len(whole_seq.split(pep)[0])
                                pos_str = ''
                                for x in pep:
                                    pos_str += str(first_pos) + '-'
                                    first_pos += 1
                                pos_str = pos_str[:-1]

                                inf_table.append([whole_seq, pep, rev_inf, pos_str])

            inf_table.sort(key=lambda x: x[2], reverse=True)
            #print(inf_table[0])
            #print(inf_table[-1])

            category = 'none'

            #print(len(inf_table))
            if len(inf_table) > 0:
                aff_1 = [inf_table[0]]
                if float(inf_table[0][2]) >= high_cutoff:
                    category = 'strong'
                elif med_cutoff <= float(inf_table[0][2]) < high_cutoff:
                    category = 'weak'

            if len(inf_table) > 1:
                aff_2 = [inf_table[1]]
            if len(inf_table) > 2:
                aff_3 = [inf_table[2]]

            vals = []
            vals_sel = []
            for pos in range(0,len(whole_seq)):
                pos_aff = 0.0
                for h in inf_table:
                    pos_split = h[3].split('-')
                    for p in pos_split:
                        if int(p) == pos:
                            pos_aff += float(h[2])
                vals.append(pos_aff)

                if len(inf_table) > 0:
                    pos_aff = 0.0
                    for h in aff_1:
                        pos_split = h[3].split('-')
                        for p in pos_split:
                            if int(p) == pos:
                                pos_aff += float(h[2])
                    vals_sel.append(pos_aff)

                if len(inf_table) > 1:
                    pos_aff = 0.0
                    for h in aff_2:
                        pos_split = h[3].split('-')
                        for p in pos_split:
                            if int(p) == pos:
                                pos_aff += float(h[2])
                    vals_sel.append(pos_aff)

                if len(inf_table) > 2:
                    pos_aff = 0.0
                    for h in aff_3:
                        pos_split = h[3].split('-')
                        for p in pos_split:
                            if int(p) == pos:
                                pos_aff += float(h[2])
                    vals_sel.append(pos_aff)


            aff_landscape = ''
            aff_1_landscape = ''
            aff_2_landscape = ''
            aff_3_landscape = ''
            for pos in range(0,len(whole_seq)):
                pos_aff = 0.0
                for h in inf_table:
                    pos_split = h[3].split('-')
                    for p in pos_split:
                        if int(p) == pos:
                            pos_aff += float(h[2])
                if len(vals) > 0 and max(vals) > 0.0:
                    pos_aff = pos_aff / float(max(vals))
                aff_landscape += str(pos_aff) + '--'

                if len(inf_table) > 0:
                    pos_aff = 0.0
                    for h in aff_1:
                        pos_split = h[3].split('-')
                        for p in pos_split:
                            if int(p) == pos:
                                pos_aff += float(h[2])
                    #pos_aff = pos_aff / float(max(vals_sel))
                    aff_1_landscape += str(pos_aff) + '--'

                if len(inf_table) > 1:
                    pos_aff = 0.0
                    for h in aff_2:
                        pos_split = h[3].split('-')
                        for p in pos_split:
                            if int(p) == pos:
                                pos_aff += float(h[2])
                    #pos_aff = pos_aff / float(max(vals_sel))
                    aff_2_landscape += str(pos_aff) + '--'

                if len(inf_table) > 2:
                    pos_aff = 0.0
                    for h in aff_3:
                        pos_split = h[3].split('-')
                        for p in pos_split:
                            if int(p) == pos:
                                pos_aff += float(h[2])
                    #pos_aff = pos_aff / float(max(vals_sel))
                    aff_3_landscape += str(pos_aff) + '--'

            aff_landscape = aff_landscape[:-2]
            if len(inf_table) > 0:
                aff_1_landscape = aff_1_landscape[:-2]
            if len(inf_table) > 1:
                aff_2_landscape = aff_2_landscape[:-2]
            if len(inf_table) > 2:
                aff_3_landscape = aff_3_landscape[:-2]


            intens_landscape = ''
            core_peps = []
            landscape_ind = ''
            pep_ind = ''
            for i in range(0,len(headers)):
                if 'Intensity Landscape in ' + exp == headers[i]:
                    landscape_ind = i
                if 'Unique Peptides in ' + exp == headers[i]:
                    pep_ind = i

            for d in raw:
#                if d[inds[2]] in whole_seq:
                if seq_start + d[inds[2]] + seq_end == whole_seq:
                    #intens_landscape = d[landscape_ind].split('--')
                    intens_landscape = d[landscape_ind]
                    peps_split = d[pep_ind].split(';')
                
                    #add_front = '0.0--' * fig_margin_left
                    #add_back = '--0.0' * fig_margin_right
                    add_front = '0.0--' * len(seq_start)
                    add_back = '--0.0' * len(seq_end)
                    intens_landscape = add_front + intens_landscape + add_back

                    landscape_split = intens_landscape.split('--')
                    #if len(landscape_split) != len(whole_seq):
                    #    print len(d[landscape_ind].split('--'))
                    #    print 'add_front', add_front
                    #    print 'add_back', add_back
                    #    print len(whole_seq)
                    #    print len(landscape_split)
                    #    print landscape_split
                    #    return None

                    pep_pos = []
                    for f in peps_split:
                        if f != '':
                            start_ind = len(whole_seq.split(f)[0])
                            pep_pos.append([start_ind, f])
                    pep_pos.sort()

                    #core_peps = peps_to_test
                    for f in pep_pos:
                        if f[1] not in core_peps:
                            core_peps.append(f[1])
                            #print(f[1])
            #print(intens_landscape)

            intens_split = intens_landscape.split('--')
            save_file = True
            #print(intens_split)
            if len(intens_split) == 0 or len(intens_split) == 1:
                intens_split = ['0.0']*len(whole_seq)
                save_file = False
                #print(intens_landscape)
            aff_split = aff_landscape.split('--')
            if len(inf_table) > 0:
                aff_1_split = aff_1_landscape.split('--')
            if len(inf_table) > 1:
                aff_2_split = aff_2_landscape.split('--')
            if len(inf_table) > 2:
                aff_3_split = aff_3_landscape.split('--')

            res_vals = []
            intens_vals = []
            aff_vals = []
            aff_1_vals = []
            aff_2_vals = []
            aff_3_vals = []
            aff_1_vals_inv = []
            aff_2_vals_inv = []
            aff_3_vals_inv = []
            core_res_vals = []
            core_intens_vals = []
            core_aff_vals = []

            all_vals_2 = []


#            print 'Core epitope', core_epi
            for g in range(0,len(whole_seq)):
                pos = g
                res = whole_seq[g]
                
                #print(len(whole_seq), len(intens_split))
                intens = intens_split[g]
                aff = aff_split[g]
                if len(inf_table) > 0:
                    aff_1 = aff_1_split[g]
                if len(inf_table) > 1:
                    aff_2 = aff_2_split[g]
                else:
                    aff_2 = 0.0
                if len(inf_table) > 2:
                    aff_3 = aff_3_split[g]
                else:
                    aff_3 = 0.0
                for i in range(0,len(core_poses)):
                    core_pos = core_poses[i]
                    if pos in core_pos:
                        if category != 'none':
                            pass
                            #print(pos, res+'*', intens, aff, aff_1, aff_2, aff_3)
                            #print(pos, res+'-'+str(i), intens, aff, aff_1, aff_2, aff_3)
                        res = res+'-'+str(i)
                else:
                    if category != 'none':
                        pass
                        #print(pos, res, intens, aff, aff_1, aff_2, aff_3)
                
                all_vals_2.append([pos, res, intens, aff, aff_1, aff_2, aff_3])
            
            fig_seq_str = whole_seq
            fig_seq = []
            for j in range(0,len(fig_seq_str)):
                fig_seq.append(fig_seq_str[j])

            res_vals = []
            intens_vals = []
            aff_vals = []
            aff_1_vals = []
            aff_2_vals = []
            aff_3_vals = []
            core_intens_vals = []
            core_aff_vals = []
            core_pep = [' ']*len(all_vals_2)
            net_pep_1 = [' ']*len(all_vals_2)
            net_pep_2 = [' ']*len(all_vals_2)
            net_pep_3 = [' ']*len(all_vals_2)
            whole_res_str = ''

            core_lists = []
            for h in range(0,len(core_poses)):
                core_lists.append([0.0]*len(all_vals_2))
            core_peps_1 = []
            for h in range(0,len(core_poses)):
                core_peps_1.append([' ']*len(all_vals_2))

            ind_count = 0
            for g in all_vals_2:
                res_draw = g[1]
                if len(g[1]) > 1:
                    #res_draw = res_draw[:1]
                    res_draw = res_draw.split('-')[0]
                res_vals.append(res_draw)
                whole_res_str += res_draw

                if g[2] != '':
                    intens_vals.append(float(g[2]))
                else:
                    intens_vals.append(0.0)
                aff_vals.append(float(g[3]))
                aff_1_vals.append(float(g[4]))
                if float(g[4]) != 0.0:
                    aff_1_vals_inv.append(1.0/float(g[4]))
                else:
                    aff_1_vals_inv.append(0.0)
                if float(g[4]) != 0.0:
                    net_pep_1[ind_count] = res_draw
                aff_2_vals.append(float(g[5]))
                if float(g[5]) != 0.0:
                    net_pep_2[ind_count] = res_draw
                aff_3_vals.append(float(g[6]))
                if float(g[6]) != 0.0:
                    net_pep_3[ind_count] = res_draw
                #if '*' in g[1]:
                if len(g[1]) > 1:
                    #ind = int(g[1][1])
                    inds_2 = g[1].split('-')[1:]
                    for ind in inds_2:
                        ind = int(ind)
                        core_lists[ind][ind_count] = float(g[2])
                        core_peps_1[ind][ind_count] = res_draw
                    #core_intens_vals.append(float(g[2]))
                    #core_aff_vals.append(float(g[3]))
                    #core_pep[ind_count] = res_draw
                #else:
                    #core_intens_vals.append(0.0)
                    #core_aff_vals.append(0.0)

                ind_count += 1

            net_label_1 = 'NetMHCIIpan #1'
            net_label_2 = 'NetMHCIIpan #2'
            net_label_3 = 'NetMHCIIpan #3'
            if max(aff_1_vals) >= high_cutoff:
                net_label_1 += ' (Strong Binding)'
            elif med_cutoff <= max(aff_1_vals) < high_cutoff:
                net_label_1 += ' (Weak Binding)'
            else:
                net_label_1 += ' (No Binding)'
            if max(aff_2_vals) >= high_cutoff:
                net_label_2 += ' (Strong Binding)'
            elif med_cutoff <= max(aff_2_vals) < high_cutoff:
                net_label_2 += ' (Weak Binding)'
            else:
                net_label_2 += ' (No Binding)'
            if max(aff_3_vals) >= high_cutoff:
                net_label_3 += ' (Strong Binding)'
            elif med_cutoff <= max(aff_3_vals) < high_cutoff:
                net_label_3 += ' (Weak Binding)'
            else:
                net_label_3 += ' (No Binding)'
            
            
            core_pep_strs = []
            for j in range(0,len(core_peps_1)):
                core_pep_str = ''
                for f in core_peps_1[j]:
                    core_pep_str += f
                core_pep_strs.append(core_pep_str)

            #core_pep_str = ''
            #for f in core_pep:
            #    core_pep_str += f
            net_pep_1_str = ''
            for f in net_pep_1:
                net_pep_1_str += f
            net_pep_2_str = ''
            for f in net_pep_2:
                net_pep_2_str += f
            net_pep_3_str = ''
            for f in net_pep_3:
                net_pep_3_str += f


            plt.rcParams.update({'font.size': 20})
            dpi = 300
            plt.rcParams['figure.dpi'] = dpi

            #plt.rcParams['font.family'] = 'serif'
            #plt.rcParams['font.serif'] = 'Ubuntu'
            #plt.rcParams['font.monospace'] = 'Ubuntu Mono'

            #font = font_manager.FontProperties(fname='/server/plateau-adm/doc/plateau.bcp.fu-berlin.de/python-script/monaco.ttf')
            #plt.rcParams['font.family'] = font.get_name()

            fig_width = 15
            fig_height = 6 # 6
            #fig = plt.figure(figsize=(fig_width,fig_height))
            fig, ax = plt.subplots(figsize=(fig_width,fig_height))
            #ax = fig.add_subplot(111)
            #ax2 = ax.twinx()
            
            #axes = [ax, ax.twinx(), ax.twinx()]
            axes = [ax, ax.twinx()]
            #fig.subplots_adjust(right=0.75)
            #axes[-1].spines['right'].set_position(('axes', 1.15))
            #axes[-1].set_frame_on(True)
            #axes[-1].patch.set_visible(False)
            #axes[-1].invert_yaxis()

            n = len(res_vals)
            ind = np.arange(n)
            width = 2.0
            y_pos = np.arange(len(res_vals))

            intens_color = '#aad7ef'
            intens_core_color = '#77a9c5'
            core_cols = [
                    '#77a9c5',
                    '#5e7eba',
                    '#5a53b2',
                    ]

            aff_color = '#e2853d'
            #aff_core_color = '#8c450f'
            aff_core_color = '#bc3212'
            aff_color_light = '#efd175'
            alpha_val_1 = 1.0
            alpha_val_2 = 1.0
            

            #ax.bar(y_pos, intens_vals, align="center", width=1, color=intens_color, alpha=alpha_val_1)
            #ax.bar(y_pos, core_intens_vals, align="center", width=1, color=intens_core_color, alpha=alpha_val_2)
            axes[0].bar(y_pos, intens_vals, align="center", width=1, color=intens_color, alpha=alpha_val_1)
            #axes[0].bar(y_pos, core_intens_vals, align="center", width=1, color=intens_core_color, alpha=alpha_val_2)
            for j in range(0,len(core_poses)):
                axes[0].bar(y_pos, core_lists[j], align="center", width=1, color=core_cols[j], alpha=alpha_val_2)

            plt.xticks(y_pos, res_vals)
            
            
            
            #ax2.step(y_pos+0.5, aff_1_vals, color=aff_core_color, lw=4)
            #ax2.step(y_pos+0.5, aff_2_vals, color=aff_color, lw=4)
            #ax2.step(y_pos+0.5, aff_3_vals, color=aff_color_light, lw=4)
            axes[1].step(y_pos+0.5, aff_1_vals, color=aff_core_color, lw=4)
            axes[1].step(y_pos+0.5, aff_2_vals, color=aff_color, lw=4)
            axes[1].step(y_pos+0.5, aff_3_vals, color=aff_color_light, lw=4)


            #axes[2].step(y_pos+0.5, aff_1_vals_inv, color=aff_core_color, lw=4, alpha=0)
            
            #ax.set_ylabel('Norm. Intens.')
            #ax2.set_ylabel('1-log50k(aff)')
            #ax2.set_ylabel('Ka [nM]')
            axes[0].set_ylabel('Norm. % Intens.')
            #ax2.set_ylabel('1-log50k(aff)')
            axes[1].set_ylabel('Ka [nM]')
            #axes[2].set_ylabel('Kd [nM]')


            med_color = '#bddb81'
            high_color = '#5aa031'
        #plt.text(0.1, high_cutoff+0.0005, 'Strong binding', color=high_color)
            #plt.axhline(y=high_cutoff, color=high_color, linewidth=2, linestyle='--')
        #plt.text(0.1, med_cutoff+0.0005, 'Weak binding', color=med_color)
        #plt.axhline(y=med_cutoff, color=med_color, linewidth=2, linestyle='--')
            
            axes[1].text(0.1, high_cutoff+0.0005, 'Strong binding', color=high_color, size=19)
            axes[1].text(len(y_pos)-7, high_cutoff+0.0005, 'Kd = ' + str(1.0/high_cutoff) + ' nM', color=high_color, size=19)
            axes[1].axhline(y=high_cutoff, color=high_color, linewidth=2, linestyle='--')
            axes[1].text(0.1, med_cutoff+0.0005, 'Weak binding', color=med_color, size=19)
            axes[1].text(len(y_pos)-7, med_cutoff+0.0005, 'Kd = ' + str(1.0/med_cutoff) + ' nM', color=med_color, size=19)
            axes[1].axhline(y=med_cutoff, color=med_color, linewidth=2, linestyle='--')
            #ax3 = ax.twinx()
            #ax3.spines['right'].set_position(('axes', 1.2))
            #ax3.set_frame_on(True)
            #ax3.patch.set_visible(True)
            
            #ax.set_ylim([0,max(intens_vals)+(max(intens_vals)*0.1)])
            axes[0].set_ylim([0,max(intens_vals)+(max(intens_vals)*0.1)])
            
            #ax2.set_ylim([0,max(aff_1_vals)+(max(aff_1_vals)*0.1)])
            if len(aff_1_vals) > 0:
                if max(aff_1_vals) <= high_cutoff:
                    #ax2.set_ylim([0, high_cutoff+0.0025])
                    #ax3.set_ylim([0, 1.0/(high_cutoff+0.0025)])
                    axes[1].set_ylim([0, high_cutoff+0.0025])
                    #axes[2].set_ylim([0, 1.0/(high_cutoff+0.0025)])
                else:
                    #ax2.set_ylim([0, max(aff_1_vals)+0.0025])
                    #ax3.set_ylim([0, max(aff_1_vals_inv)+0.0025])
                    axes[1].set_ylim([0, max(aff_1_vals)+0.0025])
                    #axes[2].set_ylim([0, max(aff_1_vals_inv)+0.0025])
            else:
                #ax2.set_ylim([0, high_cutoff+0.0025])
                #ax3.set_ylim([0, 1.0/(high_cutoff+0.0025)])
                axes[1].set_ylim([0, high_cutoff+0.0025])
                #axes[2].set_ylim([0, 1.0/(high_cutoff+0.0025)])
            
            #axes[-1].set_ylim(axes[-1].get_ylim()[::-1])
            #axes[-1].invert_yaxis()
            #ax2.set_ylim([0,1.1])
            
            
            if len(prot.split('_')) < 6:
                #prot_title = prot[0]
                prot_title = prot
            else:
                prot_title = ''
                prot_split = prot.split('_')
                for x in range(0,6):
                    prot_title += prot_split[x] + '_'
                prot_title = prot_title[:-1] + '...'


            # Draw peptides underneath

            #plt.gcf().canvas.draw()
            #ticks = [t for t in plt.gca().get_xticklabels()]
            #tick_pos = []
            #for i, t in enumerate(ticks):
                print("Label {}, data: {}".format(i, t.get_text()), t.get_window_extent(renderer=fig.canvas.get_renderer()))
                #tick_pos.append([i, t.get_text(), t.get_window_extent(renderer=fig.canvas.get_renderer())])
            #    pos = i
            #    res = t.get_text()
            #    bbox = t.get_window_extent(renderer=fig.canvas.get_renderer())
            #    pts = bbox.get_points()
             #   x0 = pts[0][0]
             #   y0 = pts[0][1]
             #   x1 = pts[1][0]
             #   y1 = pts[1][1]
             #   tick_pos.append([pos, res, x0, y0, x1, y1])

            #print(tick_pos[1], tick_pos[0])
            #res_spacing = float(tick_pos[1][2]) - float(tick_pos[0][2])

            
            xtickslocs = ax.get_xticks()
            ymin, _ = ax.get_ylim()
            #print('xticks pixel coordinates')
            tick_array = ax.transData.transform([(xtick, ymin) for xtick in xtickslocs])
            tick_pos = []
            for x in tick_array:
                tick_pos.append([x[0], x[1]])
                #print(tick_pos[-1])
            #res_spacing = int(tick_pos[1][0] - tick_pos[0][0])
            #res_spacing = int(float(res_spacing)/2.0)
            
            ##############
            #res_spacing = 13
            res_spacing = int(1.02*(0.0566*(float(len(whole_res_str))**2)) - (6.589*float(len(whole_res_str))) + 192.74) 
            
            #print(res_spacing)
            #res_spacing = int(float(res_spacing)*1.38)
            #res_spacing += 10
            x_axis_start = tick_pos[0][0]

            #peps_to_draw = [
            #        [core_pep_str, intens_core_color, 'PLAtEAU Core'],
            #        [net_pep_1_str, aff_core_color, net_label_1],
            #        [net_pep_2_str, aff_color, net_label_2],
            #        [net_pep_3_str, aff_color_light, net_label_3],
            #        ]
            
            peps_to_draw = []
            #print(len(core_poses))
            #print(core_pep_strs)
            #print(core_cols)
            for j in range(0,len(core_poses)):
                peps_to_draw.append([core_pep_strs[j], core_cols[j], 'PLAtEAU Core'])
            peps_to_draw.append([net_pep_1_str, aff_core_color, net_label_1])
            peps_to_draw.append([net_pep_2_str, aff_color, net_label_2])
            peps_to_draw.append([net_pep_3_str, aff_color_light, net_label_3])

            iedb_ex = ''
            iedb_cites = []
            cite_num = 1
            cites_seen = []
            for h in core_peps:
                iedb_out = iedb_lookup(h, dr_allele, iedb_path)
                #print('Looking up ' + h + ' in IEDB')
                # out = [pep, a[inds[9]], cite, iedb_num, qual, quant]
                if len(iedb_out) > 0:
                    cite = iedb_out[2]
                    qual = iedb_out[4]
                    quant = iedb_out[5]
                    iedb_ex = 'IEDBconf'
                    iedb_num = iedb_out[3].split('/')[-1]

                    if quant != '':
                        #pep_label = quant + ' (' + cite + ')'
                        pep_label = quant + ' (IEDB: ' + iedb_num + ') [' + str(cite_num) + ']'
                    else:
                        pep_label = '(' + qual + ') (IEDB: ' + iedb_num + ') [' + str(cite_num) + ']'
                    #print(h, pep_label)
                    
                    if cite not in cites_seen:
                        cites_seen.append(cite)
                        cite = '['+str(cite_num)+'] '+ cite
                        #iedb_cites.append([h, cite])
                        iedb_cites.append(cite)
                        cite_num += 1
                else:
                    pep_label = ''

                peps_to_draw.append([h, intens_color, pep_label])
            
            
            pep_imgs = []
            for pep_in in peps_to_draw:
                pep_file = draw_pep(pep_in[0], pep_in[1], pep_in[2], fig_width, res_spacing)
                pep_imgs.append([pep_file, re.sub(' ', '', pep_in[0]), pep_in[1], pep_in[2]])

            #img_file = pep_imgs[0]
            #im = Image.open(os.path.expanduser(img_file))
            #height = im.size[1]
            #im = np.array(im).astype(np.float) / 255
            #fig.figimage(im, 0, fig.bbox.ymax - height, zorder=10)

            #jump_size = 10
            #y_coord = fig.bbox.ymin
            #for pep_img in pep_imgs:
            #    im = Image.open(pep_img)
            #    height = im.size[1]
            #    im = np.array(im).astype(np.float) / 255
                #fig.figimage(im, 0, fig.bbox.ymax - height)
            #    fig.figimage(im, 0, y_coord - height, zorder=10)
            #    y_coord -= jump_size
            #    y_coord -= height
            #    print pep_img, height, y_coord
            
            plt.title(prot_title)
            out_file = os.path.expanduser('tmp/' + exp + '_' + dr_allele + '_' + category + '_' + prot_title + '_' + re.sub('_', '', whole_seq)[:20] + '_intens_aff_plot.png')
            if save_file == True:
                #plt.savefig(out_file, bbox_inches='tight', dpi=dpi)
                plt.savefig(out_file, dpi=dpi)
                
                plt.clf()
                print('Saving ' + out_file, count, '/', len(unique_wholes))

                imgs_to_open = [out_file]
                x_coors = [0]
                #base_x_offset = int(x_axis_start) + 230
                #base_x_offset = 650
                #res_width = 50 + (res_spacing)
                #res_width = 83
                #print('base x offset:', base_x_offset)
                #print('res width:', res_width)
                
                # linear equation: num_chars = 0.0117*(x_pos) - 5.8583
                # x_pos = (num_chars + 5.8583) / 0.0117
                
                label_space = 140
                all_widths = []
                labels = ['']
                pep_cols = ['']
                for x in pep_imgs:
                    if x[1] != '':
                        imgs_to_open.append(x[0])
                        start_res_num = len(whole_res_str.split(x[1])[0])
                        
                        x_res = float(len(whole_res_str))
                        m = (0.0594*(x_res**2)) - (7.069*x_res) + 265.28
                        b = (-0.0465*(x_res**2)) + (5.1862*x_res) + 631.02
                        extra = (-0.0151*(x_res**3)) + (1.8515*(x_res**2)) - (74.578*x_res) + 1048.1

                        #x_offset = int((float(start_res_num) + 21.684) / 0.0141) - 40
                        x_offset = int((m * float(start_res_num)) + b) + extra
                        x_coors.append(int(x_offset))
                        
                        #full_len = int((float(start_res_num+len(x[1])) + 21.684) / 0.0141) - 40
                        full_len = int((m * (float(start_res_num+len(x[1])))) + b) + extra
                        #full_len = x_offset + (len(x[1]) * (res_width))
                        if x[3] != '':
                            all_widths.append(full_len)
                        #print(x[1], start_res_num, x_offset)

                        labels.append(x[3])
                        pep_cols.append(x[2])


                images = list(map(Image.open, imgs_to_open))
                widths, heights = zip(*(i.size for i in images))

                top_margin = 20
                extra_sides = 200

                total_height = sum(heights)
                max_width = max(widths)

                background = (255,255,255)
                jump_size = 25
                cite_space = 140 * len(iedb_cites)
                new_im = Image.new('RGB', (max_width + (extra_sides), total_height + (top_margin*2) + (jump_size*len(images)) + cite_space), background)
                draw = ImageDraw.Draw(new_im)

                y_offset = 25
                for o in range(0,len(images)):
                    im = images[o]
                    im_width = widths[o]
                    extra_space = (max_width + extra_sides) - im_width
                    extra_margin = int(float(extra_space)/2.0)
                    if x_coors[o] == 0:
                        x_val = extra_margin + x_coors[o]
                    else:
                        x_val = x_coors[o]
                    new_im.paste(im, (x_val,y_offset))

                    # Label
                    if labels[o] != '':
                        font = ImageFont.truetype("DejaVuSans.ttf", 65)
                        if labels[o] == 'PLAtEAU Core' or labels[o] == net_label_1 or labels[o] == net_label_2 or labels[o] == net_label_3:
                            label_col = pep_cols[o]
                        else:
                            label_col = '#565656'
                        draw.text((max(all_widths)+label_space, y_offset),labels[o],label_col,font=font)
                    
                    y_offset += im.size[1]
                    y_offset += jump_size
                    #y_offset += 80

                if len(iedb_cites) > 0:
                    for cite in iedb_cites:
                        x_offset = 2000
                        draw.text((max(all_widths)+label_space-x_offset, y_offset), cite, '#565656', font=font)
                        y_offset += 100
                        y_offset += jump_size

                img_out = os.path.expanduser(out_file.split('.png')[0] + '_with_peps.png')
                if iedb_ex != '':
                    img_out = os.path.expanduser('tmp/pngs/' + iedb_ex + '_' + img_out.split('tmp/')[-1])
                else:
                    img_out = os.path.expanduser('tmp/pngs/' + img_out.split('tmp/')[-1])

                new_im.save(img_out)








def gen_evidence_from_peptides(pep_file):
    cols = ['Sequence', 'Proteins']
    raw, headers, inds = openfile(pep_file, cols)

    exp_inds = []
    for i in range(0,len(headers)):
        if 'Intensity ' in headers[i]:
            exp = headers[i].split('Intensity ')[-1]
            exp_inds.append([exp, i])

    headers_out = ['Sequence', 'Proteins', 'Raw file', 'Experiment', 'm/z', 'Retention time', 'Intensity']
    out = []

    for a in raw:
        seq = a[inds[0]]
        prots = a[inds[1]]

        for b in exp_inds:
            exp = b[0]
            exp_ind = b[1]

            if a[exp_ind] != '':
                if float(a[exp_ind]) > 0.0:
                    entry = [seq, prots, exp, exp, 'N/A', 'N/A', a[exp_ind]]
                    out.append(entry)

    out_file = pep_file.split('.txt')[0] + '_evidence.txt'
    savefile(out_file, out, headers_out)

    return out_file




def comb_sep(exp, exp_list):
    comb_data = []
    headers_out = []
    unique_cores = []

    for exp_name in exp_list:
        file_in = exp + '_' + exp_name + '_core_epitopes_final_renorm.txt'
        raw, headers, inds = openfile(file_in, ['Core Epitopes'])

        for b in headers:
            if b not in headers_out:
                headers_out.append(b)

        for a in raw:
            if len(a) > 0:
                unique_cores.append(a[0])

    unique_cores = list(set(unique_cores))
    unique_cores.sort()

    out = []

    count = 0
    for core in unique_cores:
        count += 1
        print(core, count, '/', len(unique_cores))
        entry = ['']*len(headers_out)

        for exp_name in exp_list:
            file_in = exp + '_' + exp_name + '_core_epitopes_final_renorm.txt'
            raw, headers, inds = openfile(file_in, ['Core Epitopes'])
            ent_ind = ''
            if len(headers) > 5:
                for i in range(0,len(headers_out)):
                    if headers_out[i] == headers[5]:
                        ent_ind = i
                        prev_ent_ind = i

                exp_vals = []
                for c in raw:
                    exp_vals.append(float(c[5]))
                    if len(c) > 0:
                        if c[0] == core:
                            entry[1] = c[1]
                            entry[2] = str(len(core))
                            entry[3] = c[3]
                            entry[4] = c[4]

                            entry[ent_ind] = c[5]
                            entry[ent_ind+1] = c[6]
                            entry[ent_ind+2] = c[7]

            if ent_ind != '':
                if entry[ent_ind] == '' and len(exp_vals) > 0:
                    entry[ent_ind] = str(min(exp_vals))
                    entry[ent_ind+1] = 'N.D.'
                    entry[ent_ind+2] = 'N.D.'
            else:
                ent_ind = prev_ent_ind + 3
                prev_ent_ind = ent_ind
                if entry[ent_ind] == '' and len(exp_vals) > 0:
                    entry[ent_ind] = str(min(exp_vals))
                    entry[ent_ind+1] = 'N.D.'
                    entry[ent_ind+2] = 'N.D.'

        out.append(entry)

    out_file = exp + '_quant_sep_final.txt'
    savefile(out_file, out, headers_out)



def quant_separately(exp, evidence_file, fasta_file):
    if evidence_file.split('_')[-1] == 'peptides.txt':
        gen_evidence_from_peptides(evidence_file)
        evidence_file = evidence_file.split('.txt')[0] + '_evidence.txt'

    exp = html.escape(exp).encode("ascii", "xmlcharrefreplace")
    evidence_file = html.escape(evidence_file).encode("ascii", "xmlcharrefreplace")
    fasta_file = html.escape(fasta_file).encode("ascii", "xmlcharrefreplace")
    evidence_file = os.path.expanduser(evidence_file)
    fasta_file = os.path.expanduser(fasta_file)
    min_epi_len = 13

    exp = str(exp)
    exp = re.sub("b'", '', exp)
    exp = re.sub("'", '', exp)
    fasta_file = str(fasta_file)
    fasta_file = re.sub("b'", '', fasta_file)
    fasta_file = re.sub("'", '', fasta_file)
    evidence_file = str(evidence_file)
    #pass_file = evidence_file.split('.txt')[0] + '_intens_norm.txt'
    #epi_file = pass_file.split('.txt')[0] + '_passpeps.txt'
    #core_file = epi_file.split('.txt')[0] + '_epitopes.txt'
    #epitope_final_file = exp + '_core_epitopes_final.txt'
    #renorm_file = epitope_final_file.split('.txt')[0] + '_renorm.txt'
    #imputation = 'no_imputation'

    # For testing filtered runs
    var_file = 'ktest_vars_final.txt'
    tech_rep_min = 2
    bio_rep_min = 1
    imputation = 'lowest_all'

    # For NetMHCIIpan
    length = 11
    dr_allele = 'DRB1_0101'
    high_cutoff = 1.0/50.0
    med_cutoff = 1.0/500.0

    start_all = timeit.default_timer()

    unique_exps = []
    raw, headers, inds = openfile(evidence_file, ['Raw file'])
    for a in raw:
        unique_exps.append(a[inds[0]])
    unique_exps = list(set(unique_exps))
    unique_exps.sort()
    prev_exp_name = 'XXXX'

    orig_evidence = evidence_file
    orig_fasta = fasta_file
    exp_orig = exp

    for exp_name in unique_exps:
        fasta_file = orig_fasta

        out = []
        for a in raw:
            if a[inds[0]] == exp_name:
                out.append(a)
        evidence_file = orig_evidence.split('.txt')[0] + '_' + exp_name + '.txt'
        savefile(evidence_file, out, headers)

        if '_' + prev_exp_name in exp:
            exp = exp.split('_'+prev_exp_name)[0]
        exp += '_' + exp_name
        prev_exp_name = exp_name
    
        pass_file = evidence_file.split('.txt')[0] + '_intens_norm.txt'
        epi_file = pass_file.split('.txt')[0] + '_passpeps.txt'
        core_file = epi_file.split('.txt')[0] + '_epitopes.txt'
        epitope_final_file = exp + '_core_epitopes_final.txt'
        renorm_file = epitope_final_file.split('.txt')[0] + '_renorm.txt'

        # 1. add areas under curve
        start = timeit.default_timer()
        add_areas_to_evidence(evidence_file)
        stop = timeit.default_timer()
        print('add_areas_to_evidence time: ', stop - start)  

        # 2. for each protein / condition, get list of matching peptides
        # this is used to generate the core epitopes for each condition
        start = timeit.default_timer()
        get_cond_peps(pass_file, fasta_file)
        #get_cond_peps_filt(pass_file, fasta_file, var_file, bio_rep_min, tech_rep_min)
        fasta_file = re.sub('_small', '', fasta_file)
        fasta_file = re.sub('_small', '', fasta_file)
        fasta_file = fasta_file.split('.fasta')[0] + '_small.fasta'
        fasta_file = re.sub('_small_small', '_small', fasta_file)
        stop = timeit.default_timer()
        print('get_cond_peps time: ', stop - start)  
        exp_end = exp.split('_')[-1]

        # 3. generate epitopes from passing peptides
        start = timeit.default_timer()
        gen_epitopes(epi_file, fasta_file, min_epi_len, min_step_size, min_epi_overlap)
        stop = timeit.default_timer()
        print('gen_epitopes: ', stop - start)  

        # 4. combine unique core epitopes
        # get total rel. intensity for each condition
        start = timeit.default_timer()
        #comb_epis(pass_file, core_file, exp)
        comb_epis_2(pass_file, core_file, exp, min_epi_len, fasta_file, epi_file)
        stop = timeit.default_timer()
        print('comb_epis: ', stop - start)  

        # 5. Renormalize after filtering, etc.
        start = timeit.default_timer()
        renorm(epitope_final_file, imputation, filt_check)
        dist_fig = gen_len_dist(exp, evidence_file, renorm_file)
        venn_html = []	
        filt_params = []
        final_output(exp, renorm_file, evidence_file, fasta_file, filt_check, dist_fig, venn_html, filt_params)
        stop = timeit.default_timer()
        print('cleanup: ', stop - start)  

        # 6. Remodel the landscapes to include jumps
        start = timeit.default_timer()
        #remodel_wholes(renorm_file, fasta_file, core_file)
        stop = timeit.default_timer()
        #print('remodel: ', stop - start)  

        # 7. NetMHCIIpan
        start = timeit.default_timer()
        #netmhciipan(renorm_file, length, dr_allele, core_file, fasta_file, high_cutoff, med_cutoff, expdir, net_path, epi_file, iedb_path)
        stop = timeit.default_timer()
        #print('netmhciipan: ', stop - start)  

        stop_all = timeit.default_timer()
        print('Total: ', stop_all - start_all)  

    comb_sep(exp_orig, unique_exps)


# Report
# Take original evidence file, compare to Plateau output
# Number of peptides / epitopes (from HTML output)
# Most important information to show...
# Differential epitopes?
# Biggest differences between conditions




##################
##################
##################
##################



banned_chars = ["'", ';', ':', '"', '/', '\\', '>', '<', '?']
#banned_chars = []
banned_str = ''
for y in banned_chars:
    banned_str += y + ' '
banned_files = ['.zip', ".htaccess", ".htpasswd", '.xml', '.aspx', '.css', '.swf', '.xhtml', '.rhtml', '.shtml', '.jsp', '.js', '.pl', '.php', '.cgi']
banned_file_str = ''
for y in banned_files:
    banned_file_str += y + ' '

#print(sys.argv)

# run: python plateau-1.0.py -exp=exp_name -evidence=evidence_file -fasta=fasta_file

# Fix: allow spaces in filenames


filt_check = sys.argv[1]


#min_epi_len = 11
min_step_size = 5
min_epi_overlap = 11

if filt_check == 'no_filt':
    exp = sys.argv[2]
    exp = html.escape(exp).encode("ascii", "xmlcharrefreplace")
    evidence_file = os.path.expanduser(sys.argv[3])
    evidence_file = html.escape(evidence_file).encode("ascii", "xmlcharrefreplace")
    fasta_file = os.path.expanduser(sys.argv[4])
    fasta_file = html.escape(fasta_file).encode("ascii", "xmlcharrefreplace")
    evidence_file = os.path.expanduser(sys.argv[3])
    fasta_file = os.path.expanduser(sys.argv[4])
    min_epi_len = int(sys.argv[5])

    evidence_check = check_file(evidence_file, banned_chars, banned_files)
    fasta_check = check_file(fasta_file, banned_chars, banned_files)
    name_check = check_file(exp, banned_chars, banned_files)
    if evidence_file.split('.')[-1].lower() != 'txt':
        evidence_check = False
    if fasta_file.split('.')[-1].lower() != 'fasta':
        fasta_check = False

    if evidence_check == False or fasta_check == False or name_check == False:
        to_print = []
        to_print.append('It seems like there is a problem with the uploaded files. Please check them.')
        to_print.append('Experiment name and file names cannot contain the following characters: ' + banned_str)
        to_print.append('Only .txt and .fasta files may be uploaded')
        html_output(to_print)
    else:	
        # Check that evidence file is formatted properly	
        cols = ['Raw file', 'm/z', 'Retention time', 'Intensity']
        raw, headers, inds = openfile(evidence_file, cols)
        if len(inds) < len(cols):
            to_print = []
            to_print.append('It seems like there is a problem with the formatting of your MaxQuant evidence.txt file. Please check it.')
            html_output(to_print)	
        else:

            # 1. add areas under curve
            add_areas_to_evidence(evidence_file)

            pass_file = evidence_file.split('.txt')[0] + '_intens_norm.txt'

            # 2. for each protein / condition, get list of matching peptides
            # this is used to generate the core epitopes for each condition
            get_cond_peps(pass_file, fasta_file)

            # 3. generate epitopes from passing peptides
            epi_file = pass_file.split('.txt')[0] + '_passpeps.txt'
            gen_epitopes(epi_file, fasta_file, min_epi_len, min_step_size, min_epi_overlap)

            # 4. combine unique core epitopes
            # get total rel. intensity for each condition
            core_file = epi_file.split('.txt')[0] + '_epitopes.txt'
            comb_epis(pass_file, core_file, exp)

            epitope_final_file = exp + '_core_epitopes_final.txt'
            # Renormalize after filtering, etc.
            imputation = 'no_imputation'
            renorm(epitope_final_file, imputation, filt_check)

            renorm_file = epitope_final_file.split('.txt')[0] + '_renorm.txt'
            dist_fig = gen_len_dist(exp, evidence_file, renorm_file)

            venn_html = []	
            filt_params = []
            final_output(exp, renorm_file, evidence_file, fasta_file, filt_check, dist_fig, venn_html, filt_params)

elif filt_check == 'yes_filt':
    exp = sys.argv[2]
    evidence_file = os.path.expanduser(sys.argv[3])
    fasta_file = os.path.expanduser(sys.argv[4])
    evidence_file = os.path.expanduser(sys.argv[3])
    fasta_file = os.path.expanduser(sys.argv[4])
    min_epi_len = int(sys.argv[5])
    template_file = os.path.expanduser(expdir+sys.argv[6].split('uploads/')[-1])
	
    evidence_check = check_file(evidence_file, banned_chars, banned_files)
    fasta_check = check_file(fasta_file, banned_chars, banned_files)
    name_check = check_file(exp, banned_chars, banned_files)
    template_check = check_file(template_file, banned_chars, banned_files)
    if evidence_file.split('.')[-1].lower() != 'txt':
        evidence_check = False
    if fasta_file.split('.')[-1].lower() != 'fasta':
        fasta_check = False


    if evidence_check == False or fasta_check == False or name_check == False:
        to_print = []
        to_print.append('It seems like there is a problem with the uploaded files. Please check them.')
        to_print.append('Experiment name and file names cannot contain the following characters: ' + banned_str)
        to_print.append('Only .txt and .fasta files may be uploaded')
        html_output(to_print)
    else:	
        # Check that evidence file is formatted properly	
        cols = ['Raw file', 'm/z', 'Retention time', 'Intensity']
        raw, headers, inds = openfile(evidence_file, cols)
        if len(inds) < len(cols):
            to_print = []
            to_print.append('It seems like there is a problem with the formatting of your MaxQuant evidence.txt file. Please check it.')
            html_output(to_print)	
        else:
            params = [[exp+';'+evidence_file+';'+fasta_file+';'+str(min_epi_len)]]
            get_exps(evidence_file, params, template_file)

elif filt_check == 'filt_ready':
    param_file = sys.argv[2]
    var_file = os.path.expanduser(sys.argv[3])
    bio_rep_min = int(sys.argv[4])
    tech_rep_min = int(sys.argv[5])
    imputation = sys.argv[6]
    filt_params = [bio_rep_min, tech_rep_min, imputation]

    params = openfile(os.path.expanduser(param_file), [])[1][0].split(';')
    exp = params[0]
    evidence_file = os.path.expanduser(params[1])
    fasta_file = os.path.expanduser(params[2])
    min_epi_len = int(params[3])

    # 1. add areas under curve
    add_areas_to_evidence(evidence_file)

    pass_file = evidence_file.split('.txt')[0] + '_intens_norm.txt'

    # 2. for each protein / condition, get list of matching peptides
    # this is used to generate the core epitopes for each condition

    get_cond_peps_filt(pass_file, fasta_file, var_file, bio_rep_min, tech_rep_min)
    pep_file_out = pass_file.split('.txt')[0] + '_pep_pass_table.txt'

    # 3. generate epitopes from passing peptides
    epi_file = pass_file.split('.txt')[0] + '_passpeps_filt.txt'
    gen_epitopes(epi_file, fasta_file, min_epi_len, min_step_size, min_epi_overlap)

    # 4. combine unique core epitopes
    # get total rel. intensity for each condition
    core_file = epi_file.split('.txt')[0] + '_epitopes.txt'
    comb_epis(pass_file, core_file, exp)

    epitope_final_file = exp + '_core_epitopes_final.txt'
    # Renormalize after filtering, etc.
    renorm(epitope_final_file, imputation, filt_check)

    renorm_file = epitope_final_file.split('.txt')[0] + '_renorm.txt'
    dist_fig = gen_len_dist(exp, evidence_file, renorm_file)

    # Venn diagrams of overlap between:
    #	- each condition filtered
    #	- bio reps of each condition
    #	- technical reps of each bio rep

    # Volcanoes
    ratio_cutoff = 2.0
    p_cutoff = 0.05
	
    var_raw, v_h, v_i = openfile(os.path.expanduser(var_file), [])
    unique_conds = []
    var_raw.insert(0,v_h)
    for a in var_raw:
        cond = a[1]
        if cond not in unique_conds:
            unique_conds.append(cond)

    venn_html = []
    if len(unique_conds) > 1:
        cond_pairs = []
        for cond1 in unique_conds:
            for cond2 in unique_conds:
                if cond1 != cond2:
                    if [cond1,cond2] not in cond_pairs and [cond2,cond1] not in cond_pairs:
                        cond_pairs.append([cond1,cond2])	


        for a in cond_pairs:
            #venn_html.append('<div style="width: 1200px; height: 90px;">')

            html_add = venn_2_conds(exp, a[0], a[1], var_raw, bio_rep_min, tech_rep_min, pep_file_out, evidence_file, renorm_file)
            for b in html_add:
                venn_html.append(b)

            #venn_html.append('</div>')

            html_add = volcano_2_conds(exp, a[0], a[1], var_raw, bio_rep_min, tech_rep_min, pep_file_out, evidence_file, renorm_file, ratio_cutoff, p_cutoff)
            for b in html_add:
                venn_html.append(b)

    final_output(exp, renorm_file, evidence_file, fasta_file, filt_check, dist_fig, venn_html, filt_params)

elif filt_check == 'test':

    #exp = 'r_test'
    #peptides = 'r_peptides.txt'
    #gen_evidence_from_peptides(peptides)
    #evidence_file = peptides.split('.txt')[0] + '_evidence.txt'
    #fasta_file = 'human_review_iso.fasta'

    exp = 'k_test'
    evidence_file = 'k_test.txt'
    fasta_file = 'k_fasta.fasta'
    #evidence_file = 'test_evidence.txt'
    #fasta_file = 'test_fasta.fasta'
    #evidence_file = 'a_evidence.txt'
    #fasta_file = 'a_fasta.fasta'
    exp = html.escape(exp).encode("ascii", "xmlcharrefreplace")
    evidence_file = html.escape(evidence_file).encode("ascii", "xmlcharrefreplace")
    fasta_file = html.escape(fasta_file).encode("ascii", "xmlcharrefreplace")
    evidence_file = os.path.expanduser(evidence_file)
    fasta_file = os.path.expanduser(fasta_file)
    min_epi_len = 13

    exp = str(exp)
    exp = re.sub("b'", '', exp)
    exp = re.sub("'", '', exp)
    fasta_file = str(fasta_file)
    fasta_file = re.sub("b'", '', fasta_file)
    fasta_file = re.sub("'", '', fasta_file)
    evidence_file = str(evidence_file)
    pass_file = evidence_file.split('.txt')[0] + '_intens_norm.txt'
    epi_file = pass_file.split('.txt')[0] + '_passpeps.txt'
    core_file = epi_file.split('.txt')[0] + '_epitopes.txt'
    epitope_final_file = exp + '_core_epitopes_final.txt'
    renorm_file = epitope_final_file.split('.txt')[0] + '_renorm.txt'
    #imputation = 'no_imputation'

    # For testing filtered runs
    var_file = 'ktest_vars_final.txt'
    tech_rep_min = 2
    bio_rep_min = 1
    imputation = 'lowest_all'

    # For NetMHCIIpan
    length = 11
    dr_allele = 'DRB1_0101'
    high_cutoff = 1.0/50.0
    med_cutoff = 1.0/500.0

    start_all = timeit.default_timer()
    # 1. add areas under curve
    start = timeit.default_timer()
    add_areas_to_evidence(evidence_file)
    stop = timeit.default_timer()
    print('add_areas_to_evidence time: ', stop - start)  

    # 2. for each protein / condition, get list of matching peptides
    # this is used to generate the core epitopes for each condition
    start = timeit.default_timer()
    #get_cond_peps(pass_file, fasta_file)
    get_cond_peps_filt(pass_file, fasta_file, var_file, bio_rep_min, tech_rep_min)
    fasta_file = fasta_file.split('.fasta')[0] + '_small.fasta'
    stop = timeit.default_timer()
    print('get_cond_peps time: ', stop - start)  

    # 3. generate epitopes from passing peptides
    start = timeit.default_timer()
    gen_epitopes(epi_file, fasta_file, min_epi_len, min_step_size, min_epi_overlap)
    stop = timeit.default_timer()
    print('gen_epitopes: ', stop - start)  

    # 4. combine unique core epitopes
    # get total rel. intensity for each condition
    start = timeit.default_timer()
    #comb_epis(pass_file, core_file, exp)
    comb_epis_2(pass_file, core_file, exp, min_epi_len, fasta_file, epi_file)
    stop = timeit.default_timer()
    print('comb_epis: ', stop - start)  

    # 5. Renormalize after filtering, etc.
    start = timeit.default_timer()
    renorm(epitope_final_file, imputation, filt_check)
    dist_fig = gen_len_dist(exp, evidence_file, renorm_file)
    venn_html = []	
    filt_params = []
    final_output(exp, renorm_file, evidence_file, fasta_file, filt_check, dist_fig, venn_html, filt_params)
    stop = timeit.default_timer()
    print('cleanup: ', stop - start)  

    # 6. Remodel the landscapes to include jumps
    start = timeit.default_timer()
    #remodel_wholes(renorm_file, fasta_file, core_file)
    stop = timeit.default_timer()
    #print('remodel: ', stop - start)  

    # 7. NetMHCIIpan
    start = timeit.default_timer()
    #netmhciipan(renorm_file, length, dr_allele, core_file, fasta_file, high_cutoff, med_cutoff, expdir, net_path, epi_file, iedb_path)
    stop = timeit.default_timer()
    #print('netmhciipan: ', stop - start)  


    stop_all = timeit.default_timer()
    print('Total: ', stop_all - start_all)  

elif filt_check == 'quant_separately':
    exp = 'k_sep'
    evidence_file = 'k_test.txt'
    fasta_file = 'k_fasta.fasta'

    quant_separately(exp, evidence_file, fasta_file)

elif filt_check == 'peptides_to_evidence':
    peptides = sys.argv[2]
    gen_evidence_from_peptides(peptides)


