
  
import sys
from sys import argv
import pandas as pd
import numpy as np
import re

### USAGE ###
# python organize_and_filter_genes.py output_folder file_location/filtered_insert_ratio_file.insert_ratios
# organizes inserts into which gene they are in, counts how many inserts per allele of each gene
# filters inserts according to number of inserts required in both alleles before we want to test that gene (change below)
# now also filters inserts according to the coefficient of variation among those inserts
# for mann whitney U test, want at least 5 inserts (I think!) in each allele
# prints out the number of inserts before and after filtering
# generates a file specific to this coevar and filtering number, added to the file name at the end: "_coevar#_filter#"
# for all inserts within a given allele, calculate the mean of all the log2(39/28) for each insert, "mean39_28"

number_inserts_per_allele_needed_to_test = 2.0
cutoff_gene_cv_ratio = 20.0 #coefficient of variation of the ratio of 39/28 for an insert

def parse_file(filename, sep='\t'): 
	df = pd.read_csv(filename, sep=sep)
	file_string = str(filename)
	p = re.compile('[^/]+[^\.](?=\.)') # matches everything before the period and after the last slash, ie. the identifier of the file
	m = p.search(file_string) # finds regex in file name
	file_identifier = m.group() # prints match in string format
	return df, file_identifier

def calc_mean_insert_ratios(grp):
	grp['mean39_28'] = np.nanmean(grp['39_28_log2_br_averaged_reads']) # #nanmean
	return grp

def calc_stdev_insert_ratios(grp):
	grp['stdev9_28'] = np.nanstd(grp['39_28_log2_br_averaged_reads']) #nanstd
	return grp

def check_for_singles(grp):
	number_alleles_found = len(grp['alleles'].unique())
	if number_alleles_found == 1:
		print grp
		print 'found one'
	return grp

class GeneObject(object):

	def __init__(self, gene_name):
		"""
		Each GeneObject is a gene (not specific to sp or sc). Keeps track of how many inserts per allele there are.
		"""
		self.gene = gene_name
		self.number_sp_inserts = 0.0
		self.number_sc_inserts = 0.0
		self.sc_ratios = []
		self.sp_ratios = []

if __name__ == '__main__':

	output_folder = argv[1].strip('/')
	files = argv[2:]

	file_id_dict = {}
	df_list = []

	for each_file in files:
		df, file_identifier = parse_file(each_file)
		df_list.append(df)
		file_id_dict[file_identifier] = df

	for file_identifier, each_df in file_id_dict.iteritems():

                print('...some pandas manipulations...')
                
                each_df = each_df[each_df['39_28_log2_br_averaged_reads'].notna()] #drop rows where there was no 28 or no 39 for calculating the logratio.
        
		each_df['gene'] = each_df.loc[:, 'annotation'].copy()
		each_df['gene'] = each_df['gene'].map(lambda x: x[2:])
		each_df['allele'] = each_df.loc[:, 'annotation'].copy()
		each_df['allele'] = each_df['allele'].map(lambda x: x[0:2])
		each_df.index.name = None #added this line for py2.7 from 2.6
		grouped = each_df.groupby('annotation', sort=False).apply(calc_mean_insert_ratios)
		grouped = each_df.groupby('annotation', sort=False).apply(calc_stdev_insert_ratios)
		grouped.set_index('gene', inplace=True, drop=False)
		unique_genes = list(grouped.gene.unique())
		gene_dict = {key: None for key in unique_genes}

		print('...building gene objects...')
		for each_gene in unique_genes:
			geneobj = GeneObject(each_gene)
			gene_dict[each_gene] = geneobj

		print('...populating gene objects with inserts...')
		for index, row in grouped.iterrows():
			row_gene_name = row['gene']
			row_allele = row['allele']
			row_mean = np.nanmean(row['39_br_averaged_reads']/row['28_br_averaged_reads']) #nanmean
			gene_object = gene_dict[row_gene_name]
			if row_allele == 'sc':
				gene_object.number_sc_inserts += 1.0
				gene_object.sc_ratios.append(row_mean)
			if row_allele == 'sp':
				gene_object.number_sp_inserts += 1.0
                                gene_object.sp_ratios.append(row_mean)
                        gene_dict[row_gene_name] = gene_object
                        #print(row_gene_name)

		filtered_genes = []

		too_few_ins=0
		too_high_cv=0

		print('...filtering genes by cv and #of inserts...')

		for each_gene, gene_info in gene_dict.iteritems():
			number_sc = float(gene_info.number_sc_inserts)
 			number_sp = float(gene_info.number_sp_inserts)
			cv_sc = np.absolute(np.nanstd(gene_info.sc_ratios)/np.nanmean(gene_info.sc_ratios))#nanmean,nanstd
			cv_sp = np.absolute(np.nanstd(gene_info.sp_ratios)/np.nanmean(gene_info.sc_ratios)) #nanmean,nanstd
			if number_sc >= number_inserts_per_allele_needed_to_test and number_sp >= number_inserts_per_allele_needed_to_test:
                                if cv_sc <= cutoff_gene_cv_ratio and cv_sp <=cutoff_gene_cv_ratio:
                                        filtered_genes.append(each_gene)
                                else:
                                    too_high_cv+=1
                        else:
                            too_few_ins+=1

		filtered_gene_df = grouped[grouped['gene'].isin(filtered_genes)]

		print('...done!')

                print 'before filtering'
		print len(grouped)
		print 'after filtering'
		print len(filtered_gene_df)

		print 'too few inserts: '+str(too_few_ins)
		print 'too high cv: '+str(too_high_cv)

                filtered_gene_df.to_csv(str(output_folder)+str('/')+str(file_identifier)+'_ratiocvfilt_'+str(cutoff_gene_cv_ratio)+'_'+str(number_inserts_per_allele_needed_to_test)+'.filtered_gene_inserts', sep='\t', index=False)
		
