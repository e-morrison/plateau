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

#expdir = sys.argv[7]
#resultdir = '../nobackup/results/'+(expdir.split('../../nobackup/uploads/')[-1])
expdir = './'
resultdir = './'


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
        fasta_str += description.split('|')[1] + ';'
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
    for prot in unique_prots:
        count += 1
        split1 = fasta_str.split(prot+';')
        fasta_ind = len(split1[0].split(';'))
        prot_seq = fasta_seqs[fasta_ind-1]

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
        fasta_str += description.split('|')[1] + ';'
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
	
    count = 0
    for prot in unique_prots:
        count += 1
        split1 = fasta_str.split(prot+';')
        fasta_ind = len(split1[0].split(';'))
        prot_seq = fasta_seqs[fasta_ind-1]

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
                #print(entry[-1], count, '/', len(unique_prots))
												
    file_out = file_in.split('.txt')[0] + '_passpeps_filt.txt'
    savefile(file_out, out, headers_out)

    file_out = file_in.split('.txt')[0] + '_missing_fasta.txt'
    savefile(file_out, no_fasta, ['Missing'])



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

				





def renorm(epi_file_in, imputation, filt_check):
    raw, headers, inds = openfile(epi_file_in, [])
    for i in range(0,len(headers)):
        if headers[i] == 'Whole Epitope':
            last_ind = i + 1

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
        for col in range(last_ind,len(headers)):
            for a in raw:
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
    for col in range(last_ind,len(headers)):
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
            #print(new_val, count, '/', len(headers)+last_ind)

        print(new_sum)

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
        if a[e_i[0]] not in unique_peps:
            unique_peps.append(a[e_i[0]])
        if a[e_i[1]] not in unique_prots:
            unique_prots.append(a[e_i[1]])

    to_print = []
    to_print.append('Contains <strong>' + str(len(unique_peps)) + '</strong> unique peptides from <strong>' + str(len(unique_prots)) + '</strong> unique proteins')

    fasta_file = str(fasta_file)
    fasta_file = re.sub("b'", '', fasta_file)
    fasta_file = re.sub("'", '', fasta_file)
    to_print.append('FASTA file used: ' + fasta_file.split('/')[-1])

    if filt_check == 'no_filt' or filt_check == 'test':
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
        if a[inds[1]] not in unique_prots:
            unique_prots.append(a[inds[1]])
        lens.append(float(a[inds[2]]))
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
        if a[e_i[0]] not in unique_peptides:
            unique_peptides.append(a[e_i[0]])	
            peptide_lens.append(len(a[e_i[0]]))

    unique_epitopes = []
    epitope_lens = []
    for a in epitopes:
        if a[ep_i[0]] not in unique_epitopes:
            unique_epitopes.append(a[ep_i[0]])	

            epitope = a[ep_i[0]]
            if epitope[0] == '*':
                epitope = epitope[1:]
            if epitope[-1] == '*':
                epitope = epitope[:-1]	

            epitope_lens.append(len(epitope))
 

    fig = plt.figure()
    ax = plt.subplot(111)

    if len(peptide_lens) > 0:
        max_len = max(peptide_lens)
        min_len = min(peptide_lens)
        num_bins = max_len - min_len

    if len(peptide_lens) > 0:
        n, bins, patches = plt.hist(peptide_lens, max(peptide_lens)-min(peptide_lens), facecolor='blue', alpha=0.5, label='Peptides')
        if len(epitope_lens) > 0:
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
    exp = 'testexp2_filt'
    exp = html.escape(exp).encode("ascii", "xmlcharrefreplace")
    evidence_file = 'k_test.txt'
    #evidence_file = 'test_evidence.txt'
    evidence_file = html.escape(evidence_file).encode("ascii", "xmlcharrefreplace")
    fasta_file = 'k_fasta.fasta'
    #fasta_file = 'test_fasta.fasta'
    fasta_file = html.escape(fasta_file).encode("ascii", "xmlcharrefreplace")
    evidence_file = os.path.expanduser(evidence_file)
    fasta_file = os.path.expanduser(fasta_file)
    min_epi_len = 13

    exp = str(exp)
    exp = re.sub("b'", '', exp)
    exp = re.sub("'", '', exp)
    evidence_file = str(evidence_file)
    pass_file = evidence_file.split('.txt')[0] + '_intens_norm.txt'
    epi_file = pass_file.split('.txt')[0] + '_passpeps.txt'
    core_file = epi_file.split('.txt')[0] + '_epitopes.txt'
    epitope_final_file = exp + '_core_epitopes_final.txt'
    renorm_file = epitope_final_file.split('.txt')[0] + '_renorm.txt'
    imputation = 'no_imputation'

    # For testing filtered runs
    var_file = 'ktest_vars_final.txt'
    tech_rep_min = 2
    bio_rep_min = 1
    imputation = 'lowest_all'

    start_all = timeit.default_timer()
    # 1. add areas under curve
    start = timeit.default_timer()
    #add_areas_to_evidence(evidence_file)
    stop = timeit.default_timer()
    print('add_areas_to_evidence time: ', stop - start)  

    # 2. for each protein / condition, get list of matching peptides
    # this is used to generate the core epitopes for each condition
    start = timeit.default_timer()
    #get_cond_peps(pass_file, fasta_file)
    #get_cond_peps_filt(pass_file, fasta_file, var_file, bio_rep_min, tech_rep_min)
    stop = timeit.default_timer()
    print('get_cond_peps time: ', stop - start)  

    # 3. generate epitopes from passing peptides
    start = timeit.default_timer()
    #gen_epitopes(epi_file, fasta_file, min_epi_len, min_step_size, min_epi_overlap)
    stop = timeit.default_timer()
    print('gen_epitopes: ', stop - start)  

    # 4. combine unique core epitopes
    # get total rel. intensity for each condition
    start = timeit.default_timer()
    #comb_epis(pass_file, core_file, exp)
    stop = timeit.default_timer()
    print('comb_epis: ', stop - start)  

    # Renormalize after filtering, etc.
    start = timeit.default_timer()
    renorm(epitope_final_file, imputation, filt_check)
    dist_fig = gen_len_dist(exp, evidence_file, renorm_file)
    venn_html = []	
    filt_params = []
    final_output(exp, renorm_file, evidence_file, fasta_file, filt_check, dist_fig, venn_html, filt_params)
    stop = timeit.default_timer()
    print('cleanup: ', stop - start)  

    stop_all = timeit.default_timer()
    print('Total: ', stop_all - start_all)  


