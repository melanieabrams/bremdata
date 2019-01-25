from sys import argv
from os import path

import math
import numpy as np
import matplotlib.pyplot as plt

import csv

from scipy import stats
import pandas
from patsy import dmatrices

import statsmodels.formula.api as sm
import statsmodels.stats.multitest as smm
import statsmodels.stats.anova.anova_lm as sma



##modified from Carly Weiss's RH_6-6-2016.py RHA script by Melanie Abrams###
        #current issues, for Fisher test:
                ##confirming that T0 is actually temperature 0 (i.e. 28C) not T1 (39C)
                ##no existing T1 pool
        #current issues, for Anova:
                ##is OLS the best way to regress?  seems like the way to do it is to make these regressions, and then feed the model into anova_lm..
                        ## could also use AnovaRM, but apparently repeated measures ANOVA ives you higher power especially if there is missing data.

        ##import data - I guess my plan is to just make a list of files in each of the conditions, but don't know how exactly to get from there to dataframe.
                #map_BLAT ALL of the new files, with assigned output names
                #input each of these fastq_pooled_reads to the program using get_files, filter_inserts

#Carly's usage
##
##### USAGE: python RH.py T0.fastq T1.fastq out_dir ###
##if len(argv) != 4:
##	print 'Usage: python RH.py T0.fastq_pooled_reads T1.fastq_pooled_reads out_directory'
##	exit()
##
##out_dir = argv[3]


#mine:

pools_28 = [ ]

pools_39 = [ ]

files_used = [pools_28, pools_39]


def get_files(names, sep='\t', cols=[]):
        '''parsing files, assumes JR output (number of pooled reads per insert is in 6th column), make sure to read in T0, then T1 in order!!'''
	files = []
	for fname in names:
		lines = []
		with open(fname) as f:
			f.readline()  # remove header
			for line in f:
				line_list = line.rstrip().split(sep)
				lines.append([line_list[i] for i in cols])
			files.append(lines)
	return files # a list of lists of each file by line

def filter_inserts(data):
        '''remove inserts mapping to the plasmid and to non-coding regions, return filtered_files'''
	filtered_files = []
	for each_file in data:
		filtered_file = []
		total_inserts = 0
		total_map_plasmid = 0
		number_inserts_not_genes = 0
		too_few_reads = 0
		for each_line in each_file:
			total_inserts += 1
			if each_line[1] == 'plasmid':  #ignore inserts mapping to plasmid
				total_map_plasmid += 1
				continue
			if each_line[4] == 'NC':  #ignore inserts not mapping to genes
				number_inserts_not_genes += 1
				continue
			else: 
				filtered_file.append(each_line)
		filtered_files.append(filtered_file)
	return filtered_files  #list of files, each is a list of lines (as lists) that pass the filtering

def calc_total_mapped_reads(data):
        '''calculates and returns the total number of mapped reads per file (including reads mapping to plasmid, noncoding etc)'''
	total_mapped_reads_both = []
	for x in data:
		total_mapped_reads = 0
		for y in x:
			total_mapped_reads = total_mapped_reads + float(y[5])
		total_mapped_reads_both.append(total_mapped_reads)
	return total_mapped_reads_both  # a list of the number of total mapped reads per file in the order they were read in (T0 first, T1 second)

def make_data_dict(filt_data):  # the main chunk of code that generates a big dictionary of everything I need
	'''data_dict = {gene name: 'gene length', [sc list[each insert in the gene is a list['insert number normalized to total reads mapped in T0',
            'same but in T1' 'rel_prop'], sp list[['', '']]}  outline of how I want things to be organized'''
	data_dict = {}
	T0_test = 0
	for y in filt_data: # for each T0 and T1
		for x in y:  # for each insert
			gene_name = x[4].strip('sc').strip('sp')
			chromosome = x[1]  # pointless?
			gene_length = float(x[8]) # self explanatory
			species = x[4][0:2]  # defines the species as either sp or sc, depending on the prefix in the gene name
			ID = str(x[0])
			rel_prop = float(x[7])
			filler_number_T0 = 1 / total_mapped_reads_both[0]  #if no reads were found for that insert in either T0 or T1, need a "filler" number to use to rep 1 read so that not dividing over 0
			filler_number_T1 = 1 / total_mapped_reads_both[1]
			if gene_name not in data_dict:
				#print 'adding', gene_name
				data_dict[gene_name] = [gene_length, []] 
			if T0_test == 0:  #parsing T0 data first
				normalized_read_count = float(x[5]) / total_mapped_reads_both[0]
				data_dict[gene_name][1].append([normalized_read_count, filler_number_T1, rel_prop, species, ID])  # set the T1 value default to 1 (normalized), replaced if a real T1 value comes along
			if T0_test == 1:  # parsing T1 data
				normalized_read_count = float(x[5]) / total_mapped_reads_both[1]				
				presence_test = 0
				for i in data_dict[gene_name][1]:  # if the insert is already in the dictionary... just add the T1 value into the list
					if i[2] == rel_prop and i[3] == species: # doublecheck that it's the exact same insert
						i[1] = normalized_read_count
						presence_test += 1						
				if presence_test == 0:
					data_dict[gene_name][1].append([filler_number_T0, normalized_read_count, rel_prop, species, ID])  # otherwise, add a new list for that insert, for inserts that appear in T1 but aren't in T0
		T0_test = T0_test + 1
	#print data_dict
	return data_dict

def how_many_inserts_per_gene(data_dictionary):
	wf = open(out_dir+"number_inserts_per_gene.txt", "w")
	wf.writelines("Gene name\tSp inserts\tSc inserts\tTotal inserts\tGene Length\n")
	
	for each_gene, values in data_dictionary.iteritems():
		
		gene_name = each_gene
		gene_length = values[0]
		total_number_inserts = len(values[1])
		total_sp_inserts = 0
		total_sc_inserts = 0
		for each_insert in values[1]:
			species = each_insert[3]
			if species == 'sp':
				total_sp_inserts += 1
			if species == 'sc':
				total_sc_inserts += 1
		wf.writelines(str(gene_name)+"\t"+str(total_sp_inserts)+"\t"+str(total_sc_inserts)+"\t"+str(total_number_inserts)+"\t"+str(gene_length)+"\n")

	wf.close()

def filter_data_dict(data_dictionary, threshold=1, lower_limit=0.0, higher_limit=1.0, number_inserts_per_allele_needed=10):
        '''filters inserts for number of reads per insert, location in gene '''
	filtered_data_dict = {}
	removed_inserts_position = {} # maybe one day I will script it so that the ones that I don't use are added to another dict. not now though ><
	T0_threshold = threshold / total_mapped_reads_both[0]  #number of normalized reads per insert required to keep that insert
	T1_threshold = threshold / total_mapped_reads_both[1]
	for key, value in data_dictionary.iteritems(): #iterate over each gene and filter inserts by read count
		gene_length = value[0]
		filtered_data_dict[key] = [gene_length, []]
		for each_insert in value[1]: # iterate over each insert
			if each_insert[0] > T0_threshold or each_insert[1] > T1_threshold:  #either the number of reads in T0 or T1 needs to pass the threshold (but it's ok if the other doesnt)
				if each_insert[2] > lower_limit and each_insert[2] < higher_limit: # filters for position within the gene, if we eventually want to do that
					log_ratio = math.log(each_insert[1]/each_insert[0], 2)  #takes the log2 ratio of T1 over T0 normalized reads, if insert passes threshold
					new_thing = each_insert
					new_thing.append(log_ratio)  # adding the log ratio to the end of the list
					filtered_data_dict[key][1].append(new_thing)  #add filtered inserts to new dictionary, not changing the old dictionary in place
	new_filtered_data_dict = {key: value for key, value in filtered_data_dict.iteritems() if len(value[1])>0}  ## make new dictionary, getting rid of keys which didn't have any inserts that pass the threshod they are just empty lists
	final_data_dict = {}
	for key, value in new_filtered_data_dict.iteritems():  # filter out genes that don't have enough inserts in each gene
		sp_inserts = 0
		sc_inserts = 0
		for each_insert in value[1]:
			if each_insert[3] == 'sp':
				sp_inserts += 1
			if each_insert[3] == 'sc':
				sc_inserts += 1
		if sp_inserts >= number_inserts_per_allele_needed and sc_inserts >= number_inserts_per_allele_needed:  #need both the number of sp and sc inserts to be over the threshold, not just one or the other
			final_data_dict[key] = value
	#interesting = final_data_dict['YKR054C']  for interogating individual genes
	#for each_insert in interesting[1]:
		#print each_insert[0], "\t", each_insert[1], "\t", each_insert[2], "\t", each_insert[3], "\t", each_insert[4], "\t", each_insert[5]
	return final_data_dict

class HemizygoteComparison(object):
        '''now doing the RH test'''
##WILL NEED TO MODIFY THIS BECAUSE OUTPUT PROBS NOT MW STAT###
	FIELDNAMES = ["Gene", "Mann Whitney Stat", "pvalue", 
				  "padjusted", "Number Sc Inserts", "Number Sp Inserts",
				  "Average Sc Fitness", "Average Sp Fitness"]

	def __init__(self, gene_name, MW_stat, pval, n_sp_inserts, n_sc_inserts, avg_sp_fitness, avg_sc_fitness):
		self.gene_name = gene_name
		self.MW_stat = MW_stat
		self.pvalue = pval
		self.adjusted_pval = None
		self.number_sp_inserts = n_sp_inserts
		self.number_sc_inserts = n_sc_inserts
		self.average_sp_fitness = avg_sp_fitness
		self.average_sc_fitness = avg_sc_fitness

	def as_dict(self):
		return {'Gene': self.gene_name, 
				'Fvalue': self.Fval,
				'pvalue': self.pvalue,
				'padjusted': self.adjusted_pval,
				'Number Sc Inserts': self.number_sc_inserts,
				'Number Sp Inserts': self.number_sp_inserts,
				'Average Sc Fitness': self.average_sc_fitness,
				'Average Sp Fitness': self.average_sp_fitness,
				}


def reciprocal_hemizygote_test(ratios_dict):
        '''executes RH test on each gene'''
	RHlist = []
	for key, value in ratios_dict.iteritems():
		gene_name = key
		sc_list = [] # list of T1/T0 ratios for inserts
		sp_list = []
		for each_insert in value[1]:  #separating p and c insert ratios into different lists
			#print each_insert
			if each_insert[3] == 'sp':
				sp_list.append(each_insert[5])
			if each_insert[3] == 'sc':
				sc_list.append(each_insert[5])
		#print sc_list, sp_list

##         Here, Carly did a Mann Whitney rank test
##		RH_test = stats.mannwhitneyu(sc_list, sp_list) # do Mann Whitney rank test on sc vs sp insert ratios, spits out a list of MW stat and pvalue
##		MW_stat = RH_test[0]
##		pvalue = RH_test[1]

                ## need to figure out how to get the data into the linear regression
                df = 'somehow import the data as a pandas dataframe
		y,x = dmatrices('Species ~ Temperature', data=df, return_type='dataframe')
		model = sm.OLS(y, x)
		#res = model.fit()
		#print(res.summary()
		
                
		RH_test = sma(model,typ=2)
		pvalue=RH_test['PR(>F)'][2] # need to confirm the number
		Fval=RH_test['F'][2]
                print RH_test
                break
		
		number_sp_inserts = len(sp_list)
		number_sc_inserts = len(sc_list)
		average_sp_fitness = sum(sp_list) / float(len(sp_list))
		average_sc_fitness = sum(sc_list) / float(len(sc_list))
		comparison_data = HemizygoteComparison(gene_name, Fval, pvalue,
			number_sp_inserts, number_sc_inserts,
			average_sp_fitness, average_sc_fitness)
		RHlist.append(comparison_data)
##	
##	pvalue_list = [cd.pvalue for cd in RHlist]
##
##	#stats2 = importr('stats') # for using rpy2
##	#adjusted_pvalues = stats2.p_adjust(FloatVector(pvalue_list), method = 'BH')  #do multiple testing correction with BH method
##	fdrbh_output = smm.multipletests(pvalue_list, method='fdr_bh')
##	adjusted_pvalues = fdrbh_output[1].tolist()
##	for cd, padjusted in zip(RHlist, adjusted_pvalues):
##		cd.adjusted_pval = padjusted
##
##
##	with open(path.join(out_dir, "testRHresults.txt"), 'w') as wf:
##		writer = csv.DictWriter(wf, HemizygoteComparison.FIELDNAMES, dialect='excel-tab')
##		writer.writeheader()
##		for cd in RHlist:
##			writer.writerow(cd.as_dict())
##
##
##


####  Running the script  ####
#files_used = argv[1:3] ##from Carly's usage
#files = get_files(files_used, cols=range(0, 9))

all_data_dict = []
for category in files_used: ##Need to figure out when this needs to split according to condition/replicate b
        files=get_files(category, cols=range(0,9))
        filtered_files = filter_inserts(files)
        total_mapped_reads_both = calc_total_mapped_reads(files)
        data_dict = make_data_dict(filtered_files)
        final_data_dict = filter_data_dict(data_dict)
        all_data_dict.append(final_data_dict)
        how_many_inserts_per_gene(final_data_dict) ## calculate how many inserts per gene in each allele that pass all the thresholds

#reciprocal_hemizygote_test(final_data_dict)
#add some function here to build dataframe or bundle that into RHT?
reciprocal_hemizygote_test(all_data_dict)





