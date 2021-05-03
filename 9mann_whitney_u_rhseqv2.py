import pandas as pd
import numpy as np
import sys
from sys import argv
import re
from scipy import stats
import statsmodels.stats.multitest as smm
import math


### USAGE ###
# python 9mann_whitney_u_uncollapsedMBA.py output_folder file.filtered_gene_inserts 
# finally, computes the Mann-Whitney U test, compares the log2(exptl/ctrl) values of inserts in the Scer allele vs the Spar allele 
# test is two-tailed (want to see if scer > spar or spar > scer), so multiply pvals by 2
# applies multiple testing correction via Benjamini Hochberg
# saves a file with each gene that was tested, pval and adjusted pval
# prints out number of genes tested



def parse_file(filename, sep='\t'): 
	df = pd.read_csv(filename, sep=sep)
	file_string = str(filename)
	p = re.compile('[^/]+[^\.](?=\.)') # matches everything before the period and after the last slash, ie. the identifier of the file
	m = p.search(file_string) # finds regex in file name
	file_identifier = m.group() # prints match in string format
	return df, file_identifier

class GeneObject(object):

	def __init__(self, df):
		self.df = df
		self.mwu_test = None
		self.two_tailed_pval = None


if __name__ == '__main__':

	output_folder = argv[1].strip('/')
	files = argv[2:]

	file_id_dict = {}
	df_list = []

	for each_file in files:
		df, file_identifier = parse_file(each_file)
		df_list.append(df)
		file_id_dict[file_identifier] = df

	for file_identifier, each_df in file_id_dict.items():
		unique_genes = each_df.gene.unique()
		unique_genes_dict = {key: None for key in unique_genes}
		another_gene_dict = {key: None for key in unique_genes}
		pval_dict = {key: None for key in unique_genes}
	
		for each_gene in unique_genes:
		    gene_df = each_df.loc[each_df['gene'] == each_gene]  # generate new df for each gene that is just its inserts from both alleles
		    sc_df = gene_df.loc[gene_df['allele'] == 'sc']
		    sp_df = gene_df.loc[gene_df['allele'] == 'sp']
		    appended_df = pd.concat([sc_df, sp_df])
		    another_gene_dict[each_gene] = appended_df
		    gene_object = GeneObject(gene_df)
		    sc_inserts_ratios = []
		    sp_inserts_ratios = []
		    sc_df = gene_df.loc[gene_df['allele'] == 'sc']
		    sp_df = gene_df.loc[gene_df['allele'] == 'sp']
		    sc_columns = sc_df.columns.values
		    sp_columns = sp_df.columns.values
		    sc_bioreps = [col for col in sc_columns if col.startswith('exptl_ctrl_log2') and col.endswith('_n_av')]
		    sp_bioreps = [col for col in sp_columns if col.startswith('exptl_ctrl_log2') and col.endswith('_n_av')]
		    for biorep in sc_bioreps:
			sc_inserts_ratios+=sc_df[biorep].tolist()
		    for biorep in sp_bioreps:
			sp_inserts_ratios+=sp_df[biorep].tolist()
		    sc_inserts_ratios[:] = [x for x in sc_inserts_ratios if str(x) != 'nan']
		    sp_inserts_ratios[:] = [x for x in sp_inserts_ratios if str(x) != 'nan']
		    if each_gene == "YHR034C":
			    print(each_gene)
			    print("sc_inserts_ratios")
			    for i in sc_inserts_ratios:
				    print(i)
			    print("sp_inserts_ratios")
			    print(len(sp_inserts_ratios), "len")
			    print(type(sp_inserts_ratios),"type")
			    for i in sp_inserts_ratios:
				    print(i)
 
		    mwu_test = stats.mannwhitneyu(sc_inserts_ratios, sp_inserts_ratios)   # mwu test
		    two_tailed_pval = mwu_test[1] * 2.0
		    gene_object.mwu_test = mwu_test
		    gene_object.two_tailed_pval = two_tailed_pval
		    unique_genes_dict[each_gene] = gene_object
		    pval_dict[each_gene] = two_tailed_pval
		    if sc_inserts_ratios == [] or sp_inserts_ratios ==[]:
			    del unique_genes_dict[each_gene]
			    del pval_dict[each_gene]
			    print(each_gene)

		gene_list = []
		pvalue_list = []
		for each_gene, each_pval in pval_dict.items():
			pvalue_list.append(each_pval)
			gene_list.append(each_gene)

		#multiple testing correction
		fdrbh_output = smm.multipletests(pvalue_list, method='fdr_bh') # benjamini hochberg method
		adjusted_pvalues = fdrbh_output[1].tolist()
		zipped_pvals = zip(pvalue_list, adjusted_pvalues)
		zipped_gene_info = zip(gene_list, zipped_pvals)
		unpacked_gene_info = [(gene, pval, adjusted_pval) for (gene, (pval, adjusted_pval)) in zipped_gene_info]
		pval_df = pd.DataFrame(unpacked_gene_info, columns=['gene', 'pval', 'adjusted_pval'])
		print(len(pval_df), 'number of genes tested')

		pval_df.sort_values(by='pval', inplace=True)
		pval_df.to_csv(str(output_folder)+str('/')+str(file_identifier)+'.mwu_test_results', sep='\t', index=False)

