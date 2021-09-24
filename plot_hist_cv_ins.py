import pandas as pd
from sys import argv
import sys
import numpy as np

import matplotlib.pyplot as plt
import seaborn as sns

### USAGE ###
# python plot_hist_cv.py

##PARAMS##

file='rerun_2.0_1.1_coeffvarANDreadcutoff_uncollapsedMBA_unpaired.insert_ratios'



number_inserts_per_allele_needed_to_test = 3.0
cutoff_gene_cv_ratio = 10.0 #coefficient of variation of the ratio of exptl/ctrl for an insert

def parse_file(filename, sep='\t'): 
	df = pd.read_csv(filename, sep=sep)

	#df=df.head(1000) #FOR TEST,UNCOMMENT FOR FULL
	return df

def calc_mean_insert_ratios(grp):
	grp['meanexptl_ctrl'] = np.nanmean(grp['exptl_ctrl_log2_br_averaged_reads']) 
	return grp

def calc_stdev_insert_ratios(grp):
	grp['stdev_exptl_ctrl'] = np.nanstd(grp['exptl_ctrl_log2_br_averaged_reads']) 
	return grp

def check_for_singles(grp):
	number_alleles_found = len(grp['alleles'].unique())
	if number_alleles_found == 1:
		print(grp)
		print('found one')
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

        each_df = parse_file(file)



        print('...some pandas manipulations...')
        
        each_df = each_df[each_df['exptl_ctrl_log2_br_averaged_reads'].notna()] #drop rows where there was no ctrl or no exptl for calculating the logratio.

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
                row_mean = np.nanmean(row['exptl_br_averaged_reads']/row['ctrl_br_averaged_reads']) 
                gene_object = gene_dict[row_gene_name]
                if row_allele == 'sc':
                        gene_object.number_sc_inserts += 1.0
                        gene_object.sc_ratios.append(row_mean)
                if row_allele == 'sp':
                        gene_object.number_sp_inserts += 1.0
                        gene_object.sp_ratios.append(row_mean)
                gene_dict[row_gene_name] = gene_object
		

        filtered_genes = []

        too_few_ins=0
        too_high_cv=0

        cvs_sc=[]
        cvs_sp=[]

        nums_sc=[]
        nums_sp=[]

        print('...filtering genes by cv and #of inserts...')

        for each_gene, gene_info in gene_dict.items():
                number_sc = float(gene_info.number_sc_inserts)
                number_sp = float(gene_info.number_sp_inserts)
                cv_sc = np.absolute(np.nanstd(gene_info.sc_ratios)/np.nanmean(gene_info.sc_ratios))
                cv_sp = np.absolute(np.nanstd(gene_info.sp_ratios)/np.nanmean(gene_info.sc_ratios))
                cvs_sc.append(cv_sc)
                cvs_sp.append(cv_sp)
                nums_sc.append(number_sc)
                nums_sp.append(number_sp)
                if number_sc >= number_inserts_per_allele_needed_to_test and number_sp >= number_inserts_per_allele_needed_to_test:
                        if cv_sc <= cutoff_gene_cv_ratio and cv_sp <=cutoff_gene_cv_ratio:
                                filtered_genes.append(each_gene)
                        else:
                                too_high_cv+=1
                else:
                        too_few_ins+=1

        filtered_gene_df = grouped[grouped['gene'].isin(filtered_genes)]

        print('...done!')

        print('before filtering')
        print(len(grouped))
        print('after filtering')
        print(len(filtered_gene_df))

        print('too few inserts: '+str(too_few_ins))
        print('too high cv: '+str(too_high_cv))

        #PLOT
        #cv
        fig, ax = plt.subplots(figsize=(9,6))
        sns.histplot(cvs_sp,label= 'cv, insert in paradoxus',color='#08A5CD')
        sns.histplot(cvs_sc,label= 'cv, insert in cerevisiae',color='#F1B629')
        plt.axvline(cutoff_gene_cv_ratio,color='red',ls='--',lw=0.5)
        #ax.set_xscale('log')
        plt.legend()
        plt.xlabel('cv')
        plt.savefig('cv_ins_hist.png')

        fig, ax = plt.subplots(figsize=(9,6))
        sns.histplot(nums_sp,label= 'num insert in paradoxus',color='#08A5CD')
        sns.histplot(nums_sc,label= 'num insert in cerevisiae',color='#F1B629')
        plt.axvline(number_inserts_per_allele_needed_to_test,color='red',ls='--',lw=0.5)
        #ax.set_xscale('log')
        plt.legend()
        plt.xlabel('# inserts per gene')
        plt.savefig('num_ins_hist.png')




		
		

		
