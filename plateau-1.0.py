#!/usr/python

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

import matplotlib.pyplot as plt
from itertools import combinations
from scipy.stats import ttest_ind
from matplotlib.markers import TICKDOWN

import shutil
from shutil import copyfile
#import cv2
import scipy

import pandas as pd
from matplotlib.colors import LogNorm
from matplotlib_venn import venn3, venn3_circles, venn2, venn2_circles

#import pymzml


def openfile(filename, cols):
	raw = list(csv.reader(open(filename, 'rU'), delimiter='\t'))
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
	if headers_out != []:
		out.insert(0,headers_out)

	with open(filename, 'w') as a_file:
		for result in out:
			result = '\t'.join(result)
			a_file.write(result + '\n')




def add_areas_to_evidence(evidence_file):
		cols = ['Raw file', 'm/z', 'Retention time', 'Intensity']
		raw, headers, inds = openfile(evidence_file, cols)

		unique_exps = []
#	headers_out = ['Run', 'Total intensity',
			#'Total area under curve', 'Area/Intens']
		headers_out = ['Run', 'Total intensity']
		out = []

		for a in raw:
				if a[inds[0]] not in unique_exps:
						unique_exps.append(a[inds[0]])

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
				print entry

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
	#			area_tot = float(b[2])

				intens_norm = ''
				if intens != '':
						intens_norm = float(intens)/intens_tot

				entry[last_ind] = str(intens)
				entry[last_ind+1] = str(intens_norm)
				print entry[0], entry[-2], entry[-1]
				out.append(entry)

		out_file = evidence_file.split('.txt')[0] + '_intens_norm.txt'
		savefile(out_file, out, headers_out)










def get_cond_peps(file_in, fasta_file):
		cols = ['Sequence', 'Proteins', 'Experiment', 'Intensity (normalized)']
		raw, headers, inds = openfile(file_in, cols)

		unique_prots = []
		unique_conds = []
		unique_bio_reps = ['']
		count = 0
		prot_seqs = []

		for a in raw:
				prots = a[inds[1]].split(';')
				count += 1
				for b in prots:
						prot = b.split('CON__')[-1]
						prot_seq = ''
						if prot not in unique_prots:
								fasta_test = False
#								print prot
								for description, sequence in fasta.read(fasta_file):
										if '|'+prot+'|' in description:
												fasta_test = True
												print 'Fasta match found:', prot, count, '/', len(raw)
												prot_seq = sequence
								if fasta_test == True:
										unique_prots.append(prot)
										prot_seqs.append([prot, prot_seq])
								else:
										print 'No FASTA match:', prot

				if a[inds[2]] not in unique_conds:
						unique_conds.append(a[inds[2]])

				#if a[inds[4]] not in unique_bio_reps:
				#		unique_bio_reps.append(a[inds[4]])

		headers_out = ['Protein', 'Experiment', 'Passing Peptides', 'FASTA seq', 'Intensities (normalized)']
		out = []
		
		for prot in unique_prots:
				for cond in unique_conds:
						for bio_rep in unique_bio_reps:
								
								pep_list = []
								intens_list = []
								for a in raw:
										if prot in a[inds[1]]:
												if cond == a[inds[2]]:
														if 1==1:
																pep_list.append(a[inds[0]])
																intens_list.append(a[inds[3]])

								pep_str = ''
								for i in pep_list:
										pep_str += i + ';'
								pep_str = pep_str[:-1]

								intens_str = ''
								for i in intens_list:
										intens_str += i + ';'
								intens_str = intens_str[:-1]

								for b in prot_seqs:
										if b[0] == prot:
												prot_seq = b[1]

								entry = [prot, cond, pep_str, prot_seq, intens_str]
								out.append(entry)
								print entry
												
		file_out = file_in.split('.txt')[0] + '_passpeps.txt'
		savefile(file_out, out, headers_out)





def gen_epitopes(file_in, fasta_file, min_epi_len, min_step_size, min_epi_overlap):
		cols = ['Protein', 'Passing Peptides', 'FASTA seq', 'Intensities (normalized)']
		raw, headers, inds = openfile(file_in, cols)

		headers_out = headers[:-1]
		headers_out.append('Core epitopes')
		headers_out.append('Whole epitopes')
#		headers_out.append('PLAtEAU Spec. Counts')
		headers_out.append('PLAtEAU Intens. (norm)')
		headers_out.append('Peptides contributing to core epitopes')
		out = []

		prot_seqs = []
		for a in raw:
				prot = a[inds[0]]
				prot_seq = a[inds[2]]
				prot_seqs.append([prot, prot_seq])

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
				prot_seq = ''
				for b in prot_seqs:
						if prot == b[0]:
								prot_seq = b[1]

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
										#print prot, tot_found[-1]

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
								#print seq_list[-1]
										
						epis = []
						start_pos = []
						end_pos = []
						for i in range(1,len(seq_list)-1):
								prev_num = float(seq_list[i-1][2])
								current_num = float(seq_list[i][2])
								next_num = float(seq_list[i+1][2])

								#if current_num != 0 or next_num != 0:
										#print i, prev_num, current_num, next_num

								if current_num == 0.0 and next_num != 0.0:
										#start_pos.append(i-3) # ******
										start_pos.append(i) # ******
										#print 'Start pos:', i+1
								if i == 1 and current_num != 0.0:
										start_pos.append(i-1) 
										#print 'Start pos:', i+1

								if current_num != 0.0 and next_num == 0.0:
										end_pos.append(i+2)
										#print 'End pos:', i
										#print ''
								if i == len(seq_list)-2 and current_num != 0.0:
										end_pos.append(i+2)

						#print '**'
						#print start_pos, end_pos

						# Problem:
						# MACARPLISVYSEKGESSGKNVTLPAVFKAPIRPDIVNFVHTNLRKNNRQPYAVSELAGHQTSAESWGTGRA whole
						# contains two cores (AVSELAGHQTSAESWGTGR), (ACARPLISVYSEKGESSGK)
						# need to split whole epitopes that have two local maxima

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
										#print y, bin_str
										split_1 = bin_str.split('0')
										split_2 = []
										for z in split_1:
												if z != '':
														split_2.append(z)
										#print y, split_1, split_2										
										if len(split_2) > 1:
												if len(split_2[0]) >= min_epi_len and len(split_2[-1]) >= min_epi_len:
														found_strs.append([y, len(split_2[0]), len(split_2[-1]), bin_str])
														#print y, bin_str
														#print found_str
														#print found_strs[-1]
								if len(found_strs) > 0:
										print found_strs
										new_min = found_strs[0][0]
										bin_str = found_strs[0][-1]
										zero_list = []
										for k in bin_str.split('1'):
												if k != '':
														zero_list.append(k)
										gap_size = len(zero_list[1])
										gap_pos_bin = len(zero_list[0]) + found_strs[0][1]
										gap_pos = epi[gap_pos_bin][1]
										print 'Gap size:', gap_size, 'starting:', gap_pos_bin, gap_pos

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
												
										#print bin_str
										#print found_strs, prot
										#print ''
										#new_min = found_strs[-1][0]
										#new_min = found_strs[0][0]
										#print ''
										#print 'Old start pos:', start_pos
										#print 'Old end pos:', end_pos
										#print 'new_min', new_min
										#bin_str = found_strs[0][-1]
										#print bin_str

										#start_ind = 0
										#for h in bin_str:
										#		if h == '1':
										#				break
										#		start_ind += 1

										#new_starts = []
										#new_ends = []

										#epi_starts = [0]
										#epi_ends = []
										#for i in range(start_ind,len(epi)-found_strs[0][2]-1):
										#		prev_num = int(bin_str[i-1])
										#		current_num = int(bin_str[i])
										#		next_num = int(bin_str[i+1])
										#		pos = epi[i][1]
										#		#print pos, current_num

										#		if current_num == 0 and next_num == 1:
										#				new_ends.append(pos)
										#				new_end_pos.append(pos)
										#				#print 'new_end:', pos

										#		if current_num == 1 and next_num == 0:
										#				new_starts.append(pos)
										#				new_start_pos.append(pos)
										#				#print 'new_start:', pos

										#		if len(new_end_pos) == 0 and i == len(epi)-found_strs[0][2]-2:
										#				new_end_pos.append(pos)
												

										#new_start_pos.append(epi[0][1])
										#new_end_pos.append(new_ends[0])

										#new_start_pos.append(new_starts[0])
										#new_end_pos.append(epi[-1][1])


						if len(new_start_pos) > 0:
								#print '**', new_start_pos
								#print new_end_pos
								for f in range(0,len(new_start_pos)):
										start_pos.append(new_start_pos[f])	
										end_pos.append(new_end_pos[f])

						start_pos.sort()
						end_pos.sort()

						#if 'P04233' in prot:
						#		print start_pos
						#		print end_pos

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
										#print epi_str, x[0]
#										if epi_str in x[0]:
										if x[0] in epi_str:
												match_pep_pos.append(x)
												#print x
												
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
								#print '**', this_pep, next_pep
										#if next_pep[1] - this_pep[1] > min_step_size and 'P04233' in prot:
										step_size = next_pep[1] - this_pep[1]
										#print '&&', this_pep, next_pep, step_size, len(this_pep[0]), step_size + min_epi_overlap 
										if next_pep[1] - this_pep[1] > min_step_size and len(this_pep[0]) <= step_size + min_epi_overlap:
												print ''
												print 'Jump:', this_pep, next_pep
												jump_test = True
												if next_pep[1]-1 not in start_pos:
														new_start_pos.append(next_pep[1]) ####
												if this_pep[1]+len(this_pep[0])-1 not in end_pos:
														new_end_pos.append(this_pep[1]+len(this_pep[0]))	####
												print 'Step size:', step_size, 'pep length:', len(this_pep[0])
												#print 'Start added:', new_start_pos[-1]
												#print 'End added:', new_end_pos[-1]

												# pep_list_1 will be peptides with start index <= this_pep[1]
												# pep_list_2 will be peptides with start index > this_pep[1]
												pep_list_1_ind_end = this_pep[1]
												pep_list_1_ind_start = start_ind
												pep_list_2_ind_start = next_pep[1]
												pep_list_2_ind_end = match_pep_pos[-1][1] + len(match_pep_pos[-1][0])
												#pep_list_2_ind = next_pep[1]
												pep_list_1 = []
												pep_list_2 = []


								#print '..', match_pep_pos
								if jump_test == True:
										for j in match_pep_pos:
												if pep_list_1_ind_start <= j[1] <= pep_list_1_ind_end:
														pep_list_1.append(j[0])
												elif pep_list_2_ind_start <= j[1] <= pep_list_2_ind_end:
														pep_list_2.append(j[0])

										#print pep_list_1
										#print pep_list_2



						if jump_test == True:
								jump_count += 1
								print 'Jumps:', jump_count, '/', total_count, ' (' + str(float(jump_count)/float(total_count)*100.0)+'%)'
								

						if len(new_start_pos) > 0:
								#print '**', new_start_pos
								#print new_end_pos
								for f in range(0,len(new_start_pos)):
										start_pos.append(new_start_pos[f])	
										end_pos.append(new_end_pos[f])

						start_pos.sort()
						end_pos.sort()
						#if 'P04233' in prot and jump_test == True:
						#		print start_pos
						#		print end_pos


						#print new_start_pos, new_end_pos
						#print start_pos, end_pos

						for i in range(0,len(start_pos)):
								epi = seq_list[start_pos[i]:end_pos[i]]
								core_start = start_pos[i]
								core_end = end_pos[i]

								#found_str = ''
								founds = []
								for x in epi:
										if float(x[2]) not in founds:
												founds.append(float(x[2]))
								#		found_str += str(x[2]) + '--'
								#found_str = found_str[:-2]
								#print prot, j[1], j[2], '(',start_pos[i],')', found_str, '(',end_pos[i],')'
								#founds_for_entry.append(found_str)
								whole_epi = ''
								for l in epi:
										whole_epi += l[0]

								if jump_test == True:
										print 'whole epi:', whole_epi

								#if prot == 'P36578':
								#if whole_epi == 'MACARPLISVYSEKGESSGKNVTLPAVFKAPIRPDIVNFVHTNLRKNNRQPYAVSELAGHQTSAESWGTGRA':
								#		print prot, j[1], j[2], '(',start_pos[i],')', found_str, '(',end_pos[i],')'

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
										#print start_pos[i], pep_list_1_ind
										#print 'List 1', pep_list_1_ind_start, start_pos[i], pep_list_1_ind_end
										#print 'List 2', pep_list_2_ind_start, start_pos[i], pep_list_2_ind_end
										if list_1_jump_test == True:
												for y in range(0,len(pep_list)):
														pep = pep_list[y]
														intens = intens_list[y]
														if pep in whole_epi and pep in pep_list_1:
																pep_match.append(pep)
																intens_match.append(intens)
																#print 'List 1 pep', pep
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
										#print prot, tot_found[-1]

										#pep_pos.sort(key=lambda x: int(x[-1]))

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
												#print seq_list_1[-1]
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
														#print k
														if float(k[2]) >= num:
																core.append(k)

												#core_str = ''
												#core_start = core[0][1]
												#core_end = core[-1][1]
												#for l in core:
												#		core_str += l[0]

												if len(core) >= min_epi_len and passing_core == False:
														passing_core = True
														core_str = ''

														core_start = core[0][1]
														core_end = core[-1][1]
														for l in core:
																core_str += l[0]
														if 'P04233' in prot and jump_test == True:
																print 'Chosen core:', core_str
																print founds_for_entry[-1]												
																
												#print whole_epi, core_str

										core_str = ''
										core = seq_list[core_start-1:core_end]
										#print 'jj', core_start, core_end, seq_list[core_start:core_end]
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
														#print k
														if float(k[2]) >= num:
																core.append(k)

												#core_str = ''
												#core_start = core[0][1]
												#core_end = core[-1][1]
												#for l in core:
												#		core_str += l[0]

												if len(core) >= min_epi_len and passing_core == False:
														passing_core = True
														core_str = ''

														core_start = core[0][1]
														core_end = core[-1][1]
														for l in core:
																core_str += l[0]
														if 'P04233' in prot and jump_test == True:
																print 'Chosen core:', core_str
																print founds_for_entry[-1]	
												#print whole_epi, core_str	

										core_str = ''
										core = seq_list[core_start-1:core_end]
										#print 'jj', core_start, core_end, seq_list[core_start:core_end]
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

				#if 'P04233' in prot:
				#		#print entry
				#		print ''
				#		print core_str
				#		print whole_str

		file_out = file_in.split('.txt')[0] + '_epitopes.txt'
		savefile(file_out, out, headers_out)





def comb_epis(evidence_file, epitope_file, exp_name):
		cols = ['Protein', 'Experiment', 'Core epitopes', 'Whole epitopes',
						'PLAtEAU Intens. (norm)', 'Peptides contributing to core epitopes'
						]
		epis, epi_headers, epi_inds = openfile(epitope_file, cols)

		cols = ['Sequence', 'Proteins', 'Experiment', 'Intensity (normalized)'
						]
		raw, headers, inds = openfile(evidence_file, cols)


		# find all unique wholes
		unique_wholes = []
		for a in epis:
				wholes = a[epi_inds[3]].split(';')
                                print wholes
				if wholes[0] != '':
						for i in wholes:
								whole = i.split('-')[1]
								if whole not in unique_wholes and whole != '':
										unique_wholes.append(whole)
										#print whole

		#print len(unique_wholes)

		# match unique cores to wholes; choose longest match
		unique_cores = []
		for a in epis:
				cores = a[epi_inds[2]].split(';')
				if cores[0] != '':
						for i in cores:
								core = i.split('-')[1]
								if core not in unique_cores and '--' not in i:
										unique_cores.append(core)
										#print core

		#print len(unique_cores)

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
								print 'core:', a, 'whole:', b

				if len(whole_lens) > 0:
						longest_len = max(whole_lens)
						for i in range(0,len(whole_matches)):
								if whole_lens[i] == longest_len:
										longest_whole = whole_matches[i]
										if longest_whole not in unique_longest_wholes:
												unique_longest_wholes.append(longest_whole)
												#print longest_whole

		#print len(unique_wholes), len(unique_longest_wholes)
								

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
								#print a, b
								if a[-1] == '*' and b[-1] == '*':
										if b[:-1] in a[:-1]:
												core_matches.append(b)
												core_lens.append(len(b))
								elif a[0] == '*' and b[0] == '*':
										if b[1:] in a[1:]:
												core_matches.append(b)
												core_lens.append(len(b))

				#print a, core_matches
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
								#print core1, core2
								if core1[0] == '*':
										core1 = core1[1:]
								if core2[0] == '*':
										core2 = core2[1:]
								if core1[-1] == '*':
										core1 = core1[:-1]
								if core2[-1] == '*':
										core2 = core2[:-1]

								n_dist_1 = a.split(core1)[0]
								n_dist_2 = a.split(core2)[0]

								if n_dist_1 > n_dist_2:
										first_core = core2
										second_core = core1
								elif n_dist_1 < n_dist_2:
										first_core = core1
										second_core = core2

								new_core = first_core + second_core[-1]
								print ''
								print '******'
								print first_core, second_core, new_core
								if first_core in new_core and second_core in new_core:
										longest_core = new_core
										if longest_cores[0][0] == '*' and longest_cores[1][0] == '*':
												longest_core = '*' + longest_core
										finals.append(longest_core)
										#print '*********'
										#print ''
										#print ''
										#print '**', first_core, second_core, new_core
								else:
										# two different cores of same length
										for x in longest_cores:
												if x not in seen:
														seen.append(x)
														#unique_longest_cores.append([x, a])
														finals.append(x)
																		#print '***', x
						for longest_core in finals:
								if longest_core not in seen:
										seen.append(longest_core)
										unique_longest_cores.append([longest_core, a])
										#if longest_core == 'SKEYFSKQK':
												#print '***', longest_core
												#print unique_longest_cores[-1]



		#print len(unique_cores), len(unique_longest_cores)



		unique_conds = []

		for a in raw:
                                #print '&&&', a[inds[2]]
				if a[inds[2]] not in unique_conds:
						unique_conds.append(a[inds[2]])

		unique_conds.sort()
                unique_bio_reps = ['']
                unique_tech_reps = ['']

		headers_out = ['Core Epitope', 'Proteins', 'Core Epitope Length', 'Whole Epitope']
		for cond in unique_conds:
				for bio_rep in unique_bio_reps:
						for tech_rep in unique_tech_reps:
								#name = cond + '_' + bio_rep + '_' + tech_rep
								name = cond
                                                                headers_out.append('% Rel. Intens. in ' + name)

		out = []

#		cols = ['Protein', 'Condition', 'Bio. Rep.', 'Core epitopes', 'Whole epitopes',
#						'PLAtEAU Spec. Counts', 'Peptides contributing to core epitopes'
#						]


		used_peps = [] # [0] = peptide, [1] = rel. intensity

		for a in unique_longest_cores:
                                #print '((', a
				entry = ['']*len(headers_out)
				core = a[0]
				whole = a[1]
				start_ind = 4

				entry[0] = core
				entry[2] = str(len(core))
				entry[3] = whole

				prots = []

				for b in epis:
						cores = b[epi_inds[2]].split(';')
						wholes = b[epi_inds[3]].split(';')

						for i in range(0,len(cores)):
								if cores[0] != '':
										core_test = cores[i].split('-')[1]
								else:
										core_test = []

								if wholes[0] != '':
										whole_test = wholes[i].split('-')[1]
								else:
										whole_test = []

								if core in whole_test:
										if b[epi_inds[0]] not in prots:
												prots.append(b[epi_inds[0]])
								if core[0] == '*':
										if core[1:] in whole_test:
												if b[epi_inds[0]] not in prots:
														prots.append(b[epi_inds[0]])
								elif core[-1] == '*':
										if core[:-1] in whole_test:
												if b[epi_inds[0]] not in prots:
														prots.append(b[epi_inds[0]])

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
						bio_inds = []
						bio_passes = 0
						for bio_rep in unique_bio_reps:
								tech_inds = []
								tech_passes = 0
								for tech_rep in unique_tech_reps:
										intensity_sum = 0.0

										for b in raw:
												if b[inds[2]] == cond:

														if b[inds[0]] in whole:
																if 1==1:
																		if b[inds[3]] != '':
																				pep_pair = [b[inds[0]], b[inds[2]]]
																				if pep_pair not in used_peps:
																						intensity_sum += float(b[inds[3]])
																						used_peps.append(pep_pair)
																						#print pep_pair
																						#if 'TPLLMQALPMGALP' in core:
																						#		print 'not used:', pep_pair
										intensity_sum *= 100.0

										if intensity_sum == 0.0:
												zeroes += 1
										else:
												tech_passes += 1
										count_total += 1

										#if core == 'KVLQLINDNTATALSYG':
										#		print '***', core, cond, bio_rep, tech_rep, kd, intensity_sum
										#if core == 'GPMGNIMIDPVLGTVGFG':
										#		print '***', core, cond, bio_rep, tech_rep, kd, intensity_sum
										entry[start_ind] = str(intensity_sum)
										tech_inds.append(start_ind)
										if start_ind not in bio_inds:
												bio_inds.append(start_ind)
										if start_ind not in all_inds:
												all_inds.append(start_ind)
										start_ind += 1

								if 1 != 1:
										for i in tech_inds:
												entry[i] = '0.0'
								else:
										bio_passes += 1
										
						if 1 != 1:
								for i in bio_inds:
										entry[i] = '0.0'
				#print entry[i]
                                zeroes = 0
				count_total = 0
				for i in all_inds:
						count_total += 1
						if entry[i] == '0.0':
								zeroes += 1

				if count_total != zeroes:
						out.append(entry)
						print entry


		seen = []
		#for a in out:
		#		if a[4] not in seen:
		#				seen.append(a[4])
		#		else:
		#				print 'Problematic duplicate whole epitope:', a[4]


		file_out = exp_name + '_core_epitopes_final.txt'
		savefile(file_out, out, headers_out)

				





def renorm(epi_file_in):
		raw, headers, inds = openfile(epi_file_in, [])
		for i in range(0,len(headers)):
				if headers[i] == 'Whole Epitope':
						last_ind = i + 1

		out = []
		for a in raw:
				entry = ['']*len(headers)
				for i in range(0,len(a)):
						entry[i] = a[i]
				out.append(entry)

		for col in range(last_ind,len(headers)):
				total_sum = 0.0
				for a in raw:
						total_sum += float(a[col])

				new_sum = 0.0
				for i in range(0,len(raw)):
						val = raw[i][col]
						new_val = float(val) / total_sum * 100.0
						#if new_val != 0.0:
						#		 print val, new_val
						new_sum += new_val
						out[i][col] = str(new_val)

				#print new_sum

		out_file = epi_file_in.split('.txt')[0] + '_renorm.txt'
		savefile(out_file, out, headers)










##################
##################
##################
##################


# run: python plateau-1.0.py -exp=exp_name -evidence=evidence_file -fasta=fasta_file

exp = sys.argv[1].split('-exp=')[-1]
evidence_file = sys.argv[2].split('-evidence=')[-1]
fasta_file = sys.argv[3].split('-fasta=')[-1]

min_epi_len = 11
min_step_size = 5
min_epi_overlap = 11


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
renorm(epitope_final_file)



