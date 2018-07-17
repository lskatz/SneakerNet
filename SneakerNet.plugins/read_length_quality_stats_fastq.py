#!/usr/local/bin/python
# Khushbu Patel | 05/29/2018
#|__This script requires Python 3.4 and modules - numpy & scipy
#|__extracts the quality string and determine the length and average quality score of each read
#|__Converts the raw values for each read set into descriptive statistics
#|__Provides descriptive stats for Read Lengths and Read Qualities, number and percentage of reads below Q30 and Ambiguous base counts
#|__Outputs separate tables for different read length buckets (150bp,250bp and 300bp)
# Usage: ./read_length_quality_and_stats_fastq.py 

import numpy as np
from scipy.stats import skew,mstats
import glob
import os
import re

# ------------------------------------------ DECLARATIONS AND INITIALIZATIONS ------------------------------------------------#
quality_scores_R1 = []
quality_scores_R2 = []
average_quality = 0
read1_length = []
read2_length = []
inserts = []
insert_sizes = []
countN1 = 0
countN2 = 0
Q1_lt_30 = 0
Q2_lt_30 = 0
R1 = []
R2 = []
Q1 = []
Q2 = []
file1 = []
file2 = []

files_149 = [] #Stores paired read files
files_249 = [] #Stores paired read files
files_299 = [] #Stores paired read files

# Following lists are to store all results for 149bp bucket
N_mean_149 = ["Mean:"]
SD_149 = ["Std_Deviation:"]
Variance_149 = ["Variance"]
median_149 = ["Median"]
Q1_149 = ["1st_Quartile:"]
Q3_149 = ["3rd_Quartile:"]
lwhisker_149 = ["Lower_whisker:"]
hwhisker_149 = ["Upper_Whisker:"]
Skew_149 = ["Skewness:"]
G_mean_149 = ["Geometric_Mean:"]

qual_N_mean_149 = ["Mean:"]
qual_SD_149 = ["Std_Deviation:"]
qual_Variance_149 = ["Variance:"]
qual_median_149 = ["Median:"]
qual_Q1_149 = ["1st_Quartile:"]
qual_Q3_149 = ["3rd_Quartile:"]
qual_lwhisker_149 = ["Lower_whisker:"]
qual_hwhisker_149 = ["Upper_Whisker:"]
qual_skew_149 = ["Skewness:"]
qual_G_mean_149 = ["Geometric_Mean:"]

# Following lists are to store all results for 249bp bucket
N_mean_249 = ["Mean:"]
SD_249 = ["Std_Deviation:"]
Variance_249 = ["Variance"]
median_249 = ["Median"]
Q1_249 = ["1st_Quartile:"]
Q3_249 = ["3rd_Quartile:"]
lwhisker_249 = ["Lower_whisker:"]
hwhisker_249 = ["Upper_Whisker:"]
Skew_249 = ["Skewness:"]
G_mean_249 = ["Geometric_Mean:"]

qual_N_mean_249 = ["Mean:"]
qual_SD_249 = ["Std_Deviation:"]
qual_Variance_249 = ["Variance:"]
qual_median_249 = ["Median:"]
qual_Q1_249 = ["1st_Quartile:"]
qual_Q3_249 = ["3rd_Quartile:"]
qual_lwhisker_249 = ["Lower_whisker:"]
qual_hwhisker_249 = ["Upper_Whisker:"]
qual_skew_249 = ["Skewness:"]
qual_G_mean_249 = ["Geometric_Mean:"]

# Following lists are to store all results for 299bp bucket
N_mean_299 = ["Mean:"]
SD_299 = ["Std_Deviation:"]
Variance_299 = ["Variance"]
median_299 = ["Median"]
Q1_299 = ["1st_Quartile:"]
Q3_299 = ["3rd_Quartile:"]
lwhisker_299 = ["Lower_whisker:"]
hwhisker_299 = ["Upper_Whisker:"]
Skew_299 = ["Skewness:"]
G_mean_299 = ["Geometric_Mean:"]

qual_N_mean_299 = ["Mean:"]
qual_SD_299 = ["Std_Deviation:"]
qual_Variance_299 = ["Variance:"]
qual_median_299 = ["Median:"]
qual_Q1_299 = ["1st_Quartile:"]
qual_Q3_299 = ["3rd_Quartile:"]
qual_lwhisker_299 = ["Lower_Whisker:"]
qual_hwhisker_299 = ["Upper_Whisker:"]
qual_skew_299 = ["Skewness:"]
qual_G_mean_299 = ["Geometric_Mean:"]

total_no_reads_149 = ["Read_count:"]
total_no_reads_249 = ["Read_count:"]
total_no_reads_299 = ["Read_count:"]
qual_lt_30_149 = ["Reads_<_Q30:"]
qual_lt_30_249 = ["Reads_<_Q30:"]
qual_lt_30_299 = ["Reads_<_Q30:"]
perc_qual_lt_30_149 = ["Percentage_reads_<_Q30"]
perc_qual_lt_30_249 = ["Percentage_reads_<_Q30"]
perc_qual_lt_30_299 = ["Percentage_reads_<_Q30"]
ambi_calls_149 = ["Ambiguous_base_calls:"]
ambi_calls_249 = ["Ambiguous_base_calls:"]
ambi_calls_299 = ["Ambiguous_base_calls:"]


R_lt_149 = ["Reads_<_149:"]
R_ge_149 = ["Reads_>=_149:"]
R_lt_249 = ["Reads_<_249:"]
R_ge_249 = ["Reads_>=_249:"]
R_lt_299 = ["Reads_<_299:"]
R_ge_299 = ["Reads_>=_299:"]

r_median = 0
i_median = 0

final_perc_R1_lt_149 = ["%_Reads_<_149:"]
final_perc_R1_ge_149 = ["%_Reads_>=_149:"]

final_perc_R1_lt_249 = ["%_Reads_<_249:"]
final_perc_R1_gt_249 = ["%_Reads_>=_249:"]

final_perc_R1_lt_299 = ["%_Reads_<_299:"]
final_perc_R1_gt_299 = ["%_Reads_>=_299:"]


final_avg_quality_lt_149 = ["Average_Quality_<_149:"]
final_avg_quality_ge_149 = ["Average_Quality_>=_149:"]
final_avg_length_lt_149 = ["Average_Length_<_149"]
final_avg_length_ge_149 = ["Average_Length_>=_149"]
final_avg_quality_lt_249 = ["Average_Quality_<_249:"]
final_avg_quality_ge_249 = ["Average_Quality_>=_249:"]
final_avg_length_lt_249 = ["Average_Length_<_249"]
final_avg_length_ge_249 = ["Average_Length_>=_249"]
final_avg_quality_lt_299 = ["Average_Quality_<_299:"]
final_avg_quality_ge_299 = ["Average_Quality_>=_299:"]
final_avg_length_lt_299 = ["Average_Length_<_299"]
final_avg_length_ge_299 = ["Average_Length_>=_299"]




# ------------------------------------------ FUNCTIONS ------------------------------------------------#
# To parse fastq file
def parseFastq(fastq_infile):
	sequences = []
	qualities = []
	
	with open(fastq_infile,"r", encoding="utf8", errors='ignore') as f:
		while True:	
			f.readline()
			seq = f.readline().rstrip()		# gets sequence line
			f.readline()
			qual = f.readline().rstrip()	# gets quality line
			if len(seq) == 0:		# if seq length is 0; reached end of file so break out of the loop
				break	
			sequences.append(seq)	# append seq to sequences list
			qualities.append(qual)	# append qual to sequences list
	
	return sequences,qualities
	


# To convert ASCII  to quality scores
def phred33toQ(qual):
	return ord(qual) - 33	# ord converts char to ASCII values and returns
	


# To calculate descriptive stats
def stats(in_array):
	a = np.array(in_array)
	mean = a.mean()
	mean = round(mean)	# rounding off
	std_dev = a.std()
	std_dev = round(std_dev)	# rounding off
	variance = np.var(a)
	variance = round(variance)	# rounding off
	Q1 = np.percentile(a,25)
	Q1 = round(Q1)	# rounding off
	median = np.percentile(a,50)
	median = round(median)	# rounding off
	Q3 = np.percentile(a,75)
	Q3 = round(Q3)	# rounding off
	skewness = skew(a)
	skewness = round(skewness)	# rounding off
	geometric_mean = mstats.gmean(a)
	geometric_mean = round(geometric_mean)	# rounding off
	
	high = []
	low = []
	IQR = Q3 - Q1
	lower = Q1 - (1.5*IQR)
	upper = Q3 - (1.5*IQR)
	
	if(min(in_array) < lower):
		low_whisker = min(in_array)
	else:
		low_whisker = min(in_array)
	
	if(max(in_array) > upper):
		high_whisker = max(in_array)
	else:
		high_whisker = max(in_array)
	
	low_whisker = round(low_whisker)	# rounding off
	high_whisker = round(high_whisker)	# rounding off
	
	return mean,std_dev,variance,Q1,median,Q3,skewness,geometric_mean,low_whisker,high_whisker
	
	
	
# Ambiguous base counts
def countN(seq):
	count = 0
	for s in seq:
		count += s.count("N")
	return count

	
# quality thresholds
def Q30(qual_list):
	count_lt_30 = 0
	for x in qual_list:
		if(x >= 0 and x < 30):
			#print(x,"<","30")	# Sanity check!
			count_lt_30 += 1
		else:
			continue
	return count_lt_30
	
	
# To get average quality scores for each read1 
def qual_score(qual):
	quality_scores = []
	read_len = []
	for Q in qual:
		score = 0
		read_len.append(len(Q))
		for val in Q:
			score += phred33toQ(val)
		average_quality = (score/len(Q))	
		quality_scores.append(average_quality)	
	return read_len,quality_scores
	
	
def print_150bp():	
	print("\n\n-----Stats_for_149_bucket---------")
	print('\t','\t'.join(files_149))
	print("Read_Length_Stats:")
	print(*lwhisker_149, sep='\t')
	print(*Q1_149, sep='\t')
	print(*median_149, sep='\t')
	print(*N_mean_149, sep='\t')
	print(*G_mean_149, sep='\t')
	print(*Q3_149, sep='\t')
	print(*hwhisker_149, sep='\t')
	print(*SD_149, sep='\t')
	print(*Variance_149, sep='\t')
	print(*Skew_149, sep='\t')
	print(*total_no_reads_149, sep='\t')
	print(*R_lt_149, sep='\t')
	print(*R_ge_149, sep='\t')
	print(*final_perc_R1_lt_149, sep='\t')
	print(*final_perc_R1_ge_149, sep='\t')
	print(*final_avg_quality_lt_149, sep='\t')
	print(*final_avg_quality_ge_149, sep='\t')
	print(*final_avg_length_lt_149, sep='\t')
	print(*final_avg_length_ge_149, sep='\t')
	print("\nRead_Quality_Stats:")
	print(*qual_lwhisker_149, sep='\t')
	print(*qual_Q1_149, sep='\t')
	print(*qual_median_149, sep='\t')
	print(*qual_N_mean_149, sep='\t')
	print(*qual_G_mean_149, sep='\t')
	print(*qual_Q3_149, sep='\t')
	print(*qual_hwhisker_149, sep='\t')
	print(*qual_SD_149, sep='\t')
	print(*qual_Variance_149, sep='\t')
	print(*qual_skew_149, sep='\t')
	print(*qual_lt_30_149, sep='\t')
	print(*perc_qual_lt_30_149, sep='\t')
	print(*ambi_calls_149, sep='\t')
	
def print_250bp():	
	print("\n\n-----Stats_for_249_bucket---------")
	print('\t','\t'.join(files_249))
	print("Read_Length_Stats:")
	print(*lwhisker_249, sep='\t')
	print(*Q1_249, sep='\t')
	print(*median_249, sep='\t')
	print(*N_mean_249, sep='\t')
	print(*G_mean_249, sep='\t')
	print(*Q3_249, sep='\t')
	print(*hwhisker_249, sep='\t')
	print(*SD_249, sep='\t')
	print(*Variance_249, sep='\t')
	print(*Skew_249, sep='\t')
	print(*total_no_reads_249, sep='\t')
	print(*R_lt_249, sep='\t')
	print(*R_ge_249, sep='\t')
	print(*final_perc_R1_lt_249, sep='\t')
	print(*final_perc_R1_gt_249, sep='\t')
	print(*final_avg_quality_lt_249, sep='\t')
	print(*final_avg_quality_ge_249, sep='\t')
	print(*final_avg_length_lt_249, sep='\t')
	print(*final_avg_length_ge_249, sep='\t')
	print("\nRead_Quality_Stats:")
	print(*qual_lwhisker_249, sep='\t')
	print(*qual_Q1_249, sep='\t')
	print(*qual_median_249, sep='\t')
	print(*qual_N_mean_249, sep='\t')
	print(*qual_G_mean_249, sep='\t')
	print(*qual_Q3_249, sep='\t')
	print(*qual_hwhisker_249, sep='\t')
	print(*qual_SD_249, sep='\t')
	print(*qual_Variance_249, sep='\t')
	print(*qual_skew_249, sep='\t')
	print(*qual_lt_30_249, sep='\t')
	print(*perc_qual_lt_30_249, sep='\t')
	print(*ambi_calls_249, sep='\t')

	
def print_300bp():	
	print("\n\n-----Stats_for_299_bucket---------")
	print('\t','\t'.join(files_299))
	print("Read_Length_Stats:")
	print(*lwhisker_299, sep='\t')
	print(*Q1_299, sep='\t')
	print(*median_299, sep='\t')
	print(*N_mean_299, sep='\t')
	print(*G_mean_299, sep='\t')
	print(*Q3_299, sep='\t')
	print(*hwhisker_299, sep='\t')
	print(*SD_299, sep='\t')
	print(*Variance_299, sep='\t')
	print(*Skew_299, sep='\t')
	print(*total_no_reads_299, sep='\t')
	print(*R_lt_299, sep='\t')
	print(*R_ge_299, sep='\t')
	print(*final_perc_R1_lt_299, sep='\t')
	print(*final_perc_R1_gt_299, sep='\t')
	print(*final_avg_quality_lt_299, sep='\t')
	print(*final_avg_quality_ge_299, sep='\t')
	print(*final_avg_length_lt_299, sep='\t')
	print(*final_avg_length_ge_299, sep='\t')
	print("\nRead_Quality_Stats:")
	print(*qual_lwhisker_299, sep='\t')
	print(*qual_Q1_299, sep='\t')
	print(*qual_median_299, sep='\t')
	print(*qual_N_mean_299, sep='\t')
	print(*qual_G_mean_299, sep='\t')
	print(*qual_Q3_299, sep='\t')
	print(*qual_hwhisker_299, sep='\t')
	print(*qual_SD_299, sep='\t')
	print(*qual_Variance_299, sep='\t')
	print(*qual_skew_299, sep='\t')
	print(*qual_lt_30_299, sep='\t')
	print(*perc_qual_lt_30_299, sep='\t')
	print(*ambi_calls_299, sep='\t')

			
	
# ---------------------------------------------------- MAIN ----------------------------------------------------------------- #	


for x in os.listdir('.'):
  if re.match('.*_R1.*.fastq$|.*_1.fastq$', x):
    file1.append(x)
	
for x in os.listdir('.'):
  if re.match('.*_R2.*.*fastq$|.*_2.fastq$', x):
    file2.append(x)
	
# sorting lists for pairs to be in the same order
file1 = sorted(file1)
file2 = sorted(file2)


for f1,f2 in zip(file1,file2):
	
	
	# command line arguments
	fastq1 = f1
	fastq2 = f2

	# Parsing fastq: function call
	seqs1,quals1 = parseFastq(fastq1)	# takes in fastq file as an input from command line and passes it as an argument to parseFastq function. Returns sequences and qualities and stores in seqs & quals
	seqs2,quals2 = parseFastq(fastq2)
		
	# total number of reads
	read_count1 = len(seqs1)
	read_count2 = len(seqs2)
	
	
	# average quality scores for each read: function call
	read1_length,quality_scores_R1 = qual_score(quals1)
	read2_length,quality_scores_R2 = qual_score(quals2)
	
	# Descriptive stats for read1 length: function call (getting the median for both R1 and R2)
	mean1,stdDev1,var1,Q1_1,r_median,Q3_1,skew1,gmean1,lwhisker1,hwhisker1 = stats(read1_length)
	mean2,stdDev2,var2,Q1_2,i_median,Q3_2,skew2,gmean2,lwhisker2,hwhisker2 = stats(read2_length)
	
	
	# Result lists
	if(hwhisker1 == 151 or hwhisker1 == 152 and hwhisker2 == 151 or hwhisker2 == 152):
		files_149.extend((f1,f2))
		
		# command line arguments
		fastq1 = f1
		fastq2 = f2

		# Parsing fastq: function call
		seqs1,quals1 = parseFastq(fastq1)	# takes in fastq file as an input from command line and passes it as an argument to parseFastq function. Returns sequences and qualities and stores in seqs & quals
		seqs2,quals2 = parseFastq(fastq2)
			
		# total number of reads
		read_count1 = len(seqs1)
		read_count2 = len(seqs2)
		
		total_no_reads_149.extend((read_count1,read_count2)) # read count


		# average quality scores for each read: function call
		read1_length,quality_scores_R1 = qual_score(quals1)
		read2_length,quality_scores_R2 = qual_score(quals2)
		
		R1_lt_149 = 0
		R1_ge_149 = 0
		R2_lt_149 = 0
		R2_ge_149 = 0
		tot_len1_ge_149 = 0
		tot_len1_lt_149 = 0
		tot_len2_lt_149 = 0
		tot_len2_ge_149 = 0
		
		for x in read1_length:
			if(x < 149):
				R1_lt_149 += 1
				tot_len1_lt_149 += x
			elif(x >= 149):
				R1_ge_149 += 1
				tot_len1_ge_149 += x
		
		for x in read2_length:
			if(x < 149):
				R2_lt_149 += 1
				tot_len2_lt_149 += x
		
			elif(x >= 149):
				R2_ge_149 += 1
				tot_len2_ge_149 += x
				
				
		R_lt_149.extend((R1_lt_149,R2_lt_149))
		R_ge_149.extend((R1_ge_149,R2_ge_149))
		
		# quality threshold function call: function call
		Q1_lt_30 = Q30(quality_scores_R1)
		Q2_lt_30 = Q30(quality_scores_R2)
		qual_lt_30_149.extend((Q1_lt_30,Q2_lt_30))

		percent_reads_lt_30_R1 = Q1_lt_30/read_count1 * 100
		percent_reads_lt_30_R2 = Q2_lt_30/read_count2 * 100
		
		# rounding off 
		percent_reads_lt_30_R1 = round(percent_reads_lt_30_R1)
		percent_reads_lt_30_R2 = round(percent_reads_lt_30_R2)
	
		perc_qual_lt_30_149.extend((percent_reads_lt_30_R1,percent_reads_lt_30_R2))

		
		# Ambiguous base function call: function call
		countN1 = countN(seqs1)
		countN2 = countN(seqs2)
		ambi_calls_149.extend((countN1,countN2))
		
		# Descriptive stats for read1 length: function call
		r_mean,r_stdDev,r_var,r_Q1,r_median,r_Q3,r_skew,r_gmean,r_lwhisker,r_hwhisker = stats(read1_length)
		i_mean,i_stdDev,i_var,i_Q1,i_median,i_Q3,i_skew,i_gmean,i_lwhisker,i_hwhisker = stats(read2_length)
		
		N_mean_149.extend((r_mean,i_mean))
		SD_149.extend((r_stdDev,i_stdDev))
		Variance_149.extend((r_var,i_var))
		median_149.extend((r_median,i_median))
		Q1_149.extend((r_Q1,i_Q1))
		Q3_149.extend((r_Q3,i_Q3))
		lwhisker_149.extend((r_lwhisker,i_lwhisker))
		hwhisker_149.extend((r_hwhisker,i_hwhisker))
		Skew_149.extend((r_skew,i_skew))
		G_mean_149.extend((r_gmean,i_gmean))

		
		# Descriptive stats for Q1 quality: function call
		q_mean,q_stdDev,q_var,q_Q1,q_median,q_Q3,q_skew,q_gmean,q_lwhisker,q_hwhisker = stats(quality_scores_R1)
		s_mean,s_stdDev,s_var,s_Q1,s_median,s_Q3,s_skew,s_gmean,s_lwhisker,s_hwhisker = stats(quality_scores_R2)
		
		qual_N_mean_149.extend((q_mean,s_mean))
		qual_SD_149.extend((q_stdDev,s_stdDev))
		qual_Variance_149.extend((q_var,s_var))
		qual_median_149.extend((q_median,s_median))
		qual_Q1_149.extend((q_Q1,s_Q1))
		qual_Q3_149.extend((q_Q3,s_Q3))
		qual_lwhisker_149.extend((q_lwhisker,s_lwhisker))
		qual_hwhisker_149.extend((q_hwhisker,s_hwhisker))
		qual_skew_149.extend((q_skew,s_skew))
		qual_G_mean_149.extend((q_gmean,s_gmean))
		
		
		# Calculating percent reads above and below 149
		perc_R1_lt_149 = (R1_lt_149/read_count1) * 100
		perc_R1_ge_149 = (R1_ge_149/read_count1) * 100
		perc_R2_lt_149 = (R2_lt_149/read_count2) * 100
		perc_R2_ge_149 = (R2_ge_149/read_count2) * 100
		
		# rounding off
		perc_R1_lt_149 = round(perc_R1_lt_149)
		perc_R1_ge_149 = round(perc_R1_ge_149)
		perc_R2_lt_149 = round(perc_R2_lt_149)
		perc_R2_ge_149 = round(perc_R2_ge_149)
		
		
		final_perc_R1_lt_149.extend((perc_R1_lt_149,perc_R2_lt_149))
		final_perc_R1_ge_149.extend((perc_R1_ge_149,perc_R2_ge_149))
		
		# Average Quality score calculation 
		avg_quality_1_le_149 = 0
		avg_quality_1_gt_149 = 0
		avg_quality_2_le_149 = 0
		avg_quality_2_gt_149 = 0
		avg_length_1_le_149 = 0
		avg_length_1_gt_149 = 0
		avg_length_2_le_149 = 0
		avg_length_2_gt_149 = 0
		tot_qual1_lt_149 = 0
		tot_qual1_ge_149 = 0
		tot_qual2_lt_149 = 0
		tot_qual2_ge_149 = 0
		
		
		for l,q in zip(read1_length,quality_scores_R1):

			if(l < 149): # for lengths le 149
				tot_qual1_lt_149 += q
			elif(l >= 149):
				tot_qual1_ge_149 += q
				
		for l,q in zip(read2_length,quality_scores_R2):

			if(l < 149): # for lengths le 149
				tot_qual2_lt_149 += q
			elif(l >= 149):
				tot_qual2_ge_149 += q
				
		
		if(R1_lt_149 == 0 and R2_lt_149 == 0):	
			avg_quality_1_le_149 = 0	
			avg_quality_2_le_149 = 0
			avg_quality_1_gt_149 = tot_qual1_ge_149 / R1_ge_149	
			avg_quality_2_gt_149 = tot_qual2_ge_149 / R2_ge_149
		
		elif(R1_ge_149 == 0 and R2_ge_149 == 0):
			avg_quality_1_le_149 = tot_qual1_lt_149 / R1_lt_149 	
			avg_quality_2_le_149 = tot_qual2_lt_149 / R2_lt_149
			avg_quality_1_gt_149 = 0	
			avg_quality_2_gt_149 = 0
		
		else:
			avg_quality_1_le_149 = tot_qual1_lt_149 / R1_lt_149 	
			avg_quality_2_le_149 = tot_qual2_lt_149 / R2_lt_149
			avg_quality_1_gt_149 = tot_qual1_ge_149 / R1_ge_149	
			avg_quality_2_gt_149 = tot_qual2_ge_149 / R2_ge_149
			
	
			
				
		# rounding off
		avg_quality_1_le_149 = round(avg_quality_1_le_149)
		avg_quality_1_gt_149 = round(avg_quality_1_gt_149)
		avg_quality_2_le_149 = round(avg_quality_2_le_149)
		avg_quality_2_gt_149 = round(avg_quality_2_gt_149)
		
		final_avg_quality_lt_149.extend((avg_quality_1_le_149,avg_quality_2_le_149))
		final_avg_quality_ge_149.extend((avg_quality_1_gt_149,avg_quality_2_gt_149))
		
		# Calculating average length of reads above and below 149
		if(R1_lt_149 == 0 and R2_lt_149 == 0):	
			avg_length_1_le_149 = 0 
			avg_length_1_gt_149 = tot_len1_ge_149/R1_ge_149
			avg_length_2_le_149 = 0
			avg_length_2_gt_149 = tot_len2_ge_149/R2_ge_149
			
		elif(R1_ge_149 == 0 and R2_ge_149 == 0):
			avg_length_1_le_149 = tot_len1_lt_149/R1_lt_149 
			avg_length_1_gt_149 = 0
			avg_length_2_le_149 = tot_len2_lt_149/R2_lt_149
			avg_length_2_gt_149 = 0
			
		else:
			avg_length_1_le_149 = tot_len1_lt_149/R1_lt_149 
			avg_length_1_gt_149 = tot_len1_ge_149/R1_ge_149
			avg_length_2_le_149 = tot_len2_lt_149/R2_lt_149
			avg_length_2_gt_149 = tot_len2_ge_149/R2_ge_149
		
		# rounding off
		avg_length_1_le_149 = round(avg_length_1_le_149)
		avg_length_1_gt_149 = round(avg_length_1_gt_149)
		avg_length_2_le_149 = round(avg_length_2_le_149)
		avg_length_2_gt_149 = round(avg_length_2_gt_149)
		
		final_avg_length_lt_149.extend((avg_length_1_le_149,avg_length_2_le_149))
		final_avg_length_ge_149.extend((avg_length_1_gt_149,avg_length_2_gt_149))
		
	
				
	elif(hwhisker1 == 251 or hwhisker1 == 252 and hwhisker2 == 251 or hwhisker2 == 252 ):
		
		files_249.extend((f1,f2))
		
		# command line arguments
		fastq1 = f1
		fastq2 = f2
	
		# Parsing fastq: function call
		seqs1,quals1 = parseFastq(fastq1)	# takes in fastq file as an input from command line and passes it as an argument to parseFastq function. Returns sequences and qualities and stores in seqs & quals
		seqs2,quals2 = parseFastq(fastq2)
			
		# total number of reads
		read_count1 = len(seqs1)
		read_count2 = len(seqs2)
		total_no_reads_249.extend((read_count1,read_count2))
		
		# average quality scores for each read: function call
		read1_length,quality_scores_R1 = qual_score(quals1)
		read2_length,quality_scores_R2 = qual_score(quals2)
		
		
		R1_lt_249 = 0
		R1_ge_249 = 0
		R2_lt_249 = 0
		R2_ge_249 = 0
		tot_len1_lt_249 = 0
		tot_len1_ge_249 = 0
		tot_len2_lt_249 = 0
		tot_len2_ge_249 = 0
		
		for x in read1_length:
			if(x < 249):
				R1_lt_249 += 1
				tot_len1_lt_249 += x
			elif(x >= 249):
				R1_ge_249 += 1
				tot_len1_ge_249 += x
		
		for x in read2_length:
			if(x < 249):
				R2_lt_249 += 1
				tot_len2_lt_249 += x
		
			elif(x >= 249):
				R2_ge_249 += 1
				tot_len2_ge_249 += x
		
		
		R_lt_249.extend((R1_lt_249,R2_lt_249))
		R_ge_249.extend((R1_ge_249,R2_ge_249))
	
		# quality threshold function call: function call
		Q1_lt_30 = Q30(quality_scores_R1)
		Q2_lt_30 = Q30(quality_scores_R2)
		
		qual_lt_30_249.extend((Q1_lt_30,Q2_lt_30))

		percent_reads_lt_30_R1 = Q1_lt_30/read_count1 * 100
		percent_reads_lt_30_R2 = Q2_lt_30/read_count2 * 100
		
		# rounding off
		percent_reads_lt_30_R1 = round(percent_reads_lt_30_R1)
		percent_reads_lt_30_R2 = round(percent_reads_lt_30_R2)

		perc_qual_lt_30_249.extend((percent_reads_lt_30_R1,percent_reads_lt_30_R2)) 

		# Ambiguous base function call: function call
		countN1 = countN(seqs1)
		countN2 = countN(seqs2)
		ambi_calls_249.extend((countN1,countN2))
		
		# Descriptive stats for read1 length: function call
		r_mean,r_stdDev,r_var,r_Q1,r_median,r_Q3,r_skew,r_gmean,r_lwhisker,r_hwhisker = stats(read1_length)
		i_mean,i_stdDev,i_var,i_Q1,i_median,i_Q3,i_skew,i_gmean,i_lwhisker,i_hwhisker = stats(read2_length)
		
		N_mean_249.extend((r_mean,i_mean))
		SD_249.extend((r_stdDev,i_stdDev))
		Variance_249.extend((r_var,i_var))
		median_249.extend((r_median,i_median))
		Q1_249.extend((r_Q1,i_Q1))
		Q3_249.extend((r_Q3,i_Q3))
		lwhisker_249.extend((r_lwhisker,i_lwhisker))
		hwhisker_249.extend((r_hwhisker,i_hwhisker))
		Skew_249.extend((r_skew,i_skew))
		G_mean_249.extend((r_gmean,i_gmean))

		
		# Descriptive stats for Q1 quality: function call
		q_mean,q_stdDev,q_var,q_Q1,q_median,q_Q3,q_skew,q_gmean,q_lwhisker,q_hwhisker = stats(quality_scores_R1)
		s_mean,s_stdDev,s_var,s_Q1,s_median,s_Q3,s_skew,s_gmean,s_lwhisker,s_hwhisker = stats(quality_scores_R2)
		
		qual_N_mean_249.extend((q_mean,s_mean))
		qual_SD_249.extend((q_stdDev,s_stdDev))
		qual_Variance_249.extend((q_var,s_var))
		qual_median_249.extend((q_median,s_median))
		qual_Q1_249.extend((q_Q1,s_Q1))
		qual_Q3_249.extend((q_Q3,s_Q3))
		qual_lwhisker_249.extend((q_lwhisker,s_lwhisker))
		qual_hwhisker_249.extend((q_hwhisker,s_hwhisker))
		qual_skew_249.extend((q_skew,s_skew))
		qual_G_mean_249.extend((q_gmean,s_gmean))
		
	
		perc_R1_lt_249 = (R1_lt_249/read_count1) * 100
		perc_R1_gt_249 = (R1_ge_249/read_count1) * 100
		perc_R2_lt_249 = (R2_lt_249/read_count2) * 100
		perc_R2_gt_249 = (R2_ge_249/read_count2) * 100
		
		# rounding off
		perc_R1_lt_249 = round(perc_R1_lt_249)
		perc_R1_gt_249 = round(perc_R1_gt_249)
		perc_R2_lt_249 = round(perc_R2_lt_249)
		perc_R2_gt_249 = round(perc_R2_gt_249)
		
				
		final_perc_R1_lt_249.extend((perc_R1_lt_249,perc_R2_lt_249))
		final_perc_R1_gt_249.extend((perc_R1_gt_249,perc_R2_gt_249))

		
		# Average Quality score calculation 
		avg_quality_1_le_249 = 0
		avg_quality_1_gt_249 = 0
		avg_quality_2_le_249 = 0
		avg_quality_2_gt_249 = 0
		avg_length_1_le_249 = 0
		avg_length_1_gt_249 = 0
		avg_length_2_le_249 = 0
		avg_length_2_gt_249 = 0
		tot_qual1_lt_249 = 0
		tot_qual1_ge_249 = 0
		tot_qual2_lt_249 = 0
		tot_qual2_ge_249 = 0

		for l,q in zip(read1_length,quality_scores_R1):

			if(l < 249): # for lengths le 249
				tot_qual1_lt_249 += q
			elif(l >= 249):
				tot_qual1_ge_249 += q
				
		for l,q in zip(read2_length,quality_scores_R2):

			if(l < 249): # for lengths le 249
				tot_qual2_lt_249 += q
			elif(l >= 249):
				tot_qual2_ge_249 += q
				
		
		# Average quality per bucket
		if(R1_lt_249 == 0 and R2_lt_249 == 0):
			avg_quality_1_le_249 = 0	# if there are no reads less than 251
			avg_quality_1_gt_249 = tot_qual1_ge_249 / R1_ge_249
			avg_quality_2_le_249 = 0	# if there are no reads less than 251
			avg_quality_2_gt_249 = tot_qual2_ge_249 / R2_ge_249
			
		elif(R1_ge_249 == 0 and R2_ge_249 == 0):
			avg_quality_1_le_249 = tot_qual1_lt_249 / R1_lt_249
			avg_quality_1_gt_249 = 0
			avg_quality_2_le_249 = tot_qual2_lt_249 / R2_lt_249
			avg_quality_2_gt_249 = 0
			
		else:
			avg_quality_1_le_249 = tot_qual1_lt_249 / R1_lt_249
			avg_quality_1_gt_249 = tot_qual1_ge_249 / R1_ge_249
			avg_quality_2_le_249 = tot_qual2_lt_249 / R2_lt_249
			avg_quality_2_gt_249 = tot_qual2_ge_249 / R2_ge_249
		
		# rounding off
		avg_quality_1_le_249 = round(avg_quality_1_le_249)
		avg_quality_1_gt_249 = round(avg_quality_1_gt_249)
		avg_quality_2_le_249 = round(avg_quality_2_le_249)
		avg_quality_2_gt_249 = round(avg_quality_2_gt_249)
	
	
		final_avg_quality_lt_249.extend((avg_quality_1_le_249,avg_quality_2_le_249))
		final_avg_quality_ge_249.extend((avg_quality_1_gt_249,avg_quality_2_gt_249))
		
		if(R1_lt_249 == 0 and R2_lt_249 == 0):
			avg_length_1_le_249 = 0
			avg_length_1_gt_249 = tot_len1_ge_249 / R1_ge_249
			avg_length_2_le_249 = 0
			avg_length_2_gt_249 = tot_len2_ge_249 / R2_ge_249
		
		elif(R1_ge_249 == 0 and R2_ge_249 == 0):
			avg_length_1_le_249 = tot_len1_lt_249 / R1_lt_249
			avg_length_1_gt_249 = 0
			avg_length_2_le_249 = tot_len2_lt_249 / R2_lt_249 
			avg_length_2_gt_249 = 0
			
		else:
			avg_length_1_le_249 = tot_len1_lt_249 / R1_lt_249
			avg_length_1_gt_249 = tot_len1_ge_249 / R1_ge_249
			avg_length_2_le_249 = tot_len2_lt_249 / R2_lt_249 
			avg_length_2_gt_249 = tot_len2_ge_249 / R2_ge_249
		
		# rounding off
		avg_length_1_le_249 = round(avg_length_1_le_249)
		avg_length_1_gt_249 = round(avg_length_1_gt_249)
		avg_length_2_le_249 = round(avg_length_2_le_249)
		avg_length_2_gt_249 = round(avg_length_2_gt_249)
		
	
		final_avg_length_lt_249.extend((avg_length_1_le_249,avg_length_2_le_249))
		final_avg_length_ge_249.extend((avg_length_1_gt_249,avg_length_2_gt_249))
	
				
	else:
		files_299.extend((f1,f2))
		
		
		# command line arguments
		fastq1 = f1
		fastq2 = f2

		# Parsing fastq: function call
		seqs1,quals1 = parseFastq(fastq1)	# takes in fastq file as an input from command line and passes it as an argument to parseFastq function. Returns sequences and qualities and stores in seqs & quals
		seqs2,quals2 = parseFastq(fastq2)
			
		# total number of reads
		read_count1 = len(seqs1)
		read_count2 = len(seqs2)
		
		total_no_reads_299.extend((read_count1,read_count2))
		
		# average quality scores for each read: function call
		read1_length,quality_scores_R1 = qual_score(quals1)
		read2_length,quality_scores_R2 = qual_score(quals2)
		
		R1_lt_299 = 0
		R1_ge_299 = 0
		R2_lt_299 = 0
		R2_ge_299 = 0
		tot_len1_lt_299 = 0
		tot_len1_ge_299 = 0
		tot_len2_lt_299 = 0
		tot_len2_ge_299 = 0
		
		for x in read1_length:
			if(x < 299):
				R1_lt_299 += 1
				tot_len1_lt_299 += x
			elif(x >= 299):
				R1_ge_299 += 1
				tot_len1_ge_299 += x
		
		for x in read2_length:
			if(x < 299):
				R2_lt_299 += 1
				tot_len2_lt_299 += x
		
			elif(x >= 299):
				R2_ge_299 += 1
				tot_len2_ge_299 += x
		
		
		
		R_lt_299.extend((R1_lt_299,R2_lt_299))
		R_ge_299.extend((R1_ge_299,R2_ge_299))
		
		# quality threshold function call: function call
		Q1_lt_30 = Q30(quality_scores_R1)
		Q2_lt_30 = Q30(quality_scores_R2)
		qual_lt_30_299.extend((Q1_lt_30,Q2_lt_30))

		percent_reads_lt_30_R1 = Q1_lt_30/len(seqs1) * 100
		percent_reads_lt_30_R2 = Q2_lt_30/len(seqs2) * 100
		
		# rounding off
		percent_reads_lt_30_R1 = round(percent_reads_lt_30_R1)
		percent_reads_lt_30_R2 = round(percent_reads_lt_30_R2)
		
		perc_qual_lt_30_299.extend((percent_reads_lt_30_R1,percent_reads_lt_30_R2))

		# Ambiguous base function call: function call
		countN1 = countN(seqs1)
		countN2 = countN(seqs2)
		ambi_calls_299.extend((countN1,countN2))
		
		# Descriptive stats for read1 length: function call
		r_mean,r_stdDev,r_var,r_Q1,r_median,r_Q3,r_skew,r_gmean,r_lwhisker,r_hwhisker = stats(read1_length)
		i_mean,i_stdDev,i_var,i_Q1,i_median,i_Q3,i_skew,i_gmean,i_lwhisker,i_hwhisker = stats(read2_length)
		
		N_mean_299.extend((r_mean,i_mean))
		SD_299.extend((r_stdDev,i_stdDev))
		Variance_299.extend((r_var,i_var))
		median_299.extend((r_median,i_median))
		Q1_299.extend((r_Q1,i_Q1))
		Q3_299.extend((r_Q3,i_Q3))
		lwhisker_299.extend((r_lwhisker,i_lwhisker))
		hwhisker_299.extend((r_hwhisker,i_hwhisker))
		Skew_299.extend((r_skew,i_skew))
		G_mean_299.extend((r_gmean,i_gmean))

		
		# Descriptive stats for Q1 quality: function call
		q_mean,q_stdDev,q_var,q_Q1,q_median,q_Q3,q_skew,q_gmean,q_lwhisker,q_hwhisker = stats(quality_scores_R1)
		s_mean,s_stdDev,s_var,s_Q1,s_median,s_Q3,s_skew,s_gmean,s_lwhisker,s_hwhisker = stats(quality_scores_R2)
		
		qual_N_mean_299.extend((q_mean,s_mean))
		qual_SD_299.extend((q_stdDev,s_stdDev))
		qual_Variance_299.extend((q_var,s_var))
		qual_median_299.extend((q_median,s_median))
		qual_Q1_299.extend((q_Q1,s_Q1))
		qual_Q3_299.extend((q_Q3,s_Q3))
		qual_lwhisker_299.extend((q_lwhisker,s_lwhisker))
		qual_hwhisker_299.extend((q_hwhisker,s_hwhisker))
		qual_skew_299.extend((q_skew,s_skew))
		qual_G_mean_299.extend((q_gmean,s_gmean))
	
		
		perc_R1_lt_299 = (R1_lt_299/read_count1) * 100
		perc_R1_gt_299 = (R1_ge_299/read_count1) * 100
		perc_R2_lt_299 = (R2_lt_299/read_count2) * 100
		perc_R2_gt_299 = (R2_ge_299/read_count2) * 100
		
		# rounding off
		perc_R1_lt_299 = round(perc_R1_lt_299)
		perc_R1_gt_299 = round(perc_R1_gt_299)
		perc_R2_lt_299 = round(perc_R2_lt_299)
		perc_R2_gt_299 = round(perc_R2_gt_299)
		
		final_perc_R1_lt_299.extend((perc_R1_lt_299,perc_R2_lt_299))
		final_perc_R1_gt_299.extend((perc_R1_gt_299,perc_R2_gt_299))
		
		#header.append("\n\n-----Stats for 299 bucket---------")
		
		avg_quality_1_le_299 = 0
		avg_quality_1_gt_299 = 0
		avg_quality_2_le_299 = 0
		avg_quality_2_gt_299 = 0
		avg_length_1_le_299 = 0
		avg_length_1_gt_299 = 0
		avg_length_2_le_299 = 0
		avg_length_2_gt_299 = 0
		tot_qual1_lt_299 = 0
		tot_qual1_ge_299 = 0
		tot_qual2_lt_299 = 0
		tot_qual2_ge_299 = 0
		

		for l,q in zip(read1_length,quality_scores_R1):

			if(l <= 299): # for lengths le 249
				tot_qual1_lt_299 += q
			elif(l > 299):
				tot_qual1_ge_299 += q
				
		for l,q in zip(read2_length,quality_scores_R2):

			if(l <= 299): # for lengths le 249
				tot_qual2_lt_299 += q
			elif(l > 299):
				tot_qual2_ge_299 += q
				
	
		if(R1_lt_299 == 0 and R2_lt_299 == 0):	
			avg_quality_1_le_299 = 0
			avg_quality_1_gt_299 = tot_qual1_ge_299 / R1_ge_299
			avg_quality_2_le_299 = 0
			avg_quality_2_gt_299 = tot_qual2_ge_299 / R2_ge_299
			
		elif(R1_ge_299 == 0 and R2_ge_299 == 0):
			avg_quality_1_le_299 = tot_qual1_lt_299 / R1_lt_299
			avg_quality_1_gt_299 = 0
			avg_quality_2_le_299 = tot_qual2_lt_299 / R2_lt_299
			avg_quality_2_gt_299 = 0
			
		else:
			avg_quality_1_le_299 = tot_qual1_lt_299 / R1_lt_299
			avg_quality_1_gt_299 = tot_qual1_ge_299 / R1_ge_299
			avg_quality_2_le_299 = tot_qual2_lt_299 / R2_lt_299
			avg_quality_2_gt_299 = tot_qual2_ge_299 / R2_ge_299
		
		# rounding off upto 5 decimal places
		avg_quality_1_le_299 = round(avg_quality_1_le_299)
		avg_quality_1_gt_299 = round(avg_quality_1_gt_299)
		avg_quality_2_le_299 = round(avg_quality_2_le_299)
		avg_quality_2_gt_299 = round(avg_quality_2_gt_299)
		
		final_avg_quality_lt_299.extend((avg_quality_1_le_299,avg_quality_2_le_299))
		final_avg_quality_ge_299.extend((avg_quality_1_gt_299,avg_quality_2_gt_299))

		if(R1_lt_299 == 0 and R2_lt_299 == 0):	
			avg_length_1_le_299 = 0
			avg_length_1_gt_299 = tot_len1_ge_299 / R1_ge_299
			avg_length_2_le_299 = 0
			avg_length_2_gt_299 = tot_len2_ge_299 / R2_ge_299
			
		elif(R1_ge_299 == 0 and R2_ge_299 == 0):
			avg_length_1_le_299 = tot_len1_lt_299 / R1_lt_299
			avg_length_1_gt_299 = 0
			avg_length_2_le_299 = tot_len2_lt_299 / R2_lt_299
			avg_length_2_gt_299 = 0
			
		else:
			avg_length_1_le_299 = tot_len1_lt_299 / R1_lt_299
			avg_length_1_gt_299 = tot_len1_ge_299 / R1_ge_299
			avg_length_2_le_299 = tot_len2_lt_299 / R2_lt_299
			avg_length_2_gt_299 = tot_len2_ge_299 / R2_ge_299
		
		
		# rounding off
		avg_length_1_le_299 = round(avg_length_1_le_299)
		avg_length_1_gt_299 = round(avg_length_1_gt_299)
		avg_length_2_le_299 = round(avg_length_2_le_299)
		avg_length_2_gt_299 = round(avg_length_2_gt_299)
		
		final_avg_length_lt_299.extend((avg_length_1_le_299,avg_length_2_le_299))
		final_avg_length_ge_299.extend((avg_length_1_gt_299,avg_length_2_gt_299))

		

#function call
print_150bp()
print_250bp()
print_300bp()
