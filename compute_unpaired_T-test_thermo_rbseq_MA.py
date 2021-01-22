import sys
from sys import argv
import pandas as pd
import numpy as np
import os
from Bio import SeqIO
import operator
from datetime import datetime
import collections
import matplotlib.pyplot as plt
import argparse
import re
from scipy import stats
import statsmodels.stats.multitest as smm
import math
import matplotlib.pyplot as plt

#usage: modify parameters, run in terminal as python compute_rbseq.py 36_genefit_forTstat_MA.txt 37_genefit_fortstat_MA.txt...

#Note: modified from Jeff Skerker iPYNB notebook by mabrams

#Classes and functions
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

class GeneObject_RHpair(object):

    def __init__(self, df):
        self.df = df
        self.mwu_test = None
        self.two_tailed_pval = None

#Functions
def count_values(grp):
    return grp


def group_and_filter(biorep_fitness,temp,dropNA=True):
    grouped = biorep_fitness.groupby('annotation', sort=False).apply(count_values)
    grouped.set_index('gene', inplace=True, drop=False)
    unique_genes = list(grouped.gene.unique())
    gene_dict = {key: None for key in unique_genes}
    for each_gene in unique_genes:
        geneobj = GeneObject(each_gene)
        gene_dict[each_gene] = geneobj
    for index, row in grouped.iterrows():
        row_gene_name = row['gene']
        row_allele = row['allele']
        gene_object = gene_dict[row_gene_name]
        if row_allele == 'sc':
            gene_object.number_sc_inserts += 1.0
        if row_allele == 'sp':
            gene_object.number_sp_inserts += 1.0
            gene_dict[row_gene_name] = gene_object

    filtered_genes = []

    for each_gene, gene_info in gene_dict.items():
        number_sc = float(gene_info.number_sc_inserts)
        number_sp = float(gene_info.number_sp_inserts)
        # print(each_gene, number_sc, number_sp)
        if (number_sc == 1.0) and (number_sp == 1.0):
            filtered_genes.append(each_gene)

    filtered_gene_df = grouped[grouped['gene'].isin(filtered_genes)]

    print(len(grouped))
    print(len(filtered_gene_df))

    if dropNA==True:
        filtered_gene_df.to_csv('geneFit_'+temp+'_for_Tstat_filtered_gene_inserts.txt', sep='\t', index=False)
        filtered_gene_df.dropna(axis=0, how='any',inplace=True)
        filtered_gene_df.to_csv('geneFit_'+temp+'_dropNA.txt', sep='\t', index=False)

def calc_pairs_and_compute_stats(munged_file):
    '''computes rh pairs and stat test'''
    # first, load fitness data from munged file, then filter/dropNA

    temp=munged_file[:3]
    print(temp)
    biorep_fitness = pd.read_csv(munged_file,header=0,sep='\t',low_memory=False)
    group_and_filter(biorep_fitness,temp,dropNA=True)


    # re-load the data without the NA to create final list of RH pairs
    biorep_fitness = pd.read_csv('geneFit_'+temp+'_dropNA.txt',header=0,sep='\t',low_memory=False)
    group_and_filter(biorep_fitness,temp,dropNA=False)
    RHpair_df = pd.read_csv('geneFit_'+temp+'_dropNA.txt',header=0,sep='\t',low_memory=False)
    RHpair_df.sort_values(by=['gene','allele'],inplace=True)
    print('RHpair_df')
    print(RHpair_df)


    # RH-pair tests


    output_folder = './'
    unique_genes = RHpair_df.gene.unique()
    unique_genes_dict = {key: None for key in unique_genes}
    another_gene_dict = {key: None for key in unique_genes}
    pval_dict = {key: None for key in unique_genes}

    for each_gene in unique_genes:
        gene_df = RHpair_df.loc[RHpair_df['gene'] == each_gene]  # generate new df for each gene that is just its inserts from both alleles
        sc_df = gene_df.loc[gene_df['allele'] == 'sc']
        sp_df = gene_df.loc[gene_df['allele'] == 'sp']
        appended_df = pd.concat([sc_df, sp_df])
        another_gene_dict[each_gene] = appended_df
        gene_object = GeneObject_RHpair(gene_df)
        sc_inserts_ratios = []
        sp_inserts_ratios = []
        sc_df = gene_df.loc[gene_df['allele'] == 'sc']
        sp_df = gene_df.loc[gene_df['allele'] == 'sp']
        sc_columns = sc_df.columns.values
        sp_columns = sp_df.columns.values
        sc_bioreps = [col for col in sc_columns if col.startswith('Biorep')]
        sp_bioreps = [col for col in sp_columns if col.startswith('Biorep')]
        for biorep in sc_bioreps:
            sc_inserts_ratios+=sc_df[biorep].tolist()
        for biorep in sp_bioreps:
            sp_inserts_ratios+=sp_df[biorep].tolist()
        sc_inserts_ratios[:] = [x for x in sc_inserts_ratios if str(x) != 'nan']
        sp_inserts_ratios[:] = [x for x in sp_inserts_ratios if str(x) != 'nan']
        if each_gene == "YHR023W":
            appended_df.to_csv(str(output_folder)+str('/')+temp+'_YHR023W.txt', sep='\t',index=None)
            print(each_gene)
            print("sc_inserts_ratios")
            for i in sc_inserts_ratios:
                print(i)
            print("sp_inserts_ratios")
            print(len(sp_inserts_ratios), "len")
            print(type(sp_inserts_ratios),"type")
            for i in sp_inserts_ratios:
                print(i)
     
        
        # mwu_test = stats.wilcoxon(sc_inserts_ratios, sp_inserts_ratios)         # wilcoxon
        try:
            mwu_test = stats.mannwhitneyu(sc_inserts_ratios, sp_inserts_ratios)   # mwu test
        except ValueError:
            mwu_test = [1.0,1.0] #throws this error if they're all the same
        # mwu_test = stats.ttest_ind(sc_inserts_ratios, sp_inserts_ratios)      # t-test
        
        two_tailed_pval = mwu_test[1] * 2.0
        # two_tailed_pval = mwu_test[1]
        
        gene_object.mwu_test = mwu_test
        gene_object.two_tailed_pval = two_tailed_pval
        unique_genes_dict[each_gene] = gene_object
        pval_dict[each_gene] = two_tailed_pval
        if sc_inserts_ratios == [] or sp_inserts_ratios ==[]:
            del unique_genes_dict[each_gene]
            del pval_dict[each_gene]
            # print(each_gene)

    gene_list = []
    pvalue_list = []
    for each_gene, each_pval in pval_dict.items():
        pvalue_list.append(each_pval)
        gene_list.append(each_gene)


    print(pval_dict['YHR023W'])


    #multiple testing correction
    fdrbh_output = smm.multipletests(pvalue_list, method='fdr_bh') # benjamini hochberg method
    adjusted_pvalues = fdrbh_output[1].tolist()
    zipped_pvals = zip(pvalue_list, adjusted_pvalues)
    zipped_gene_info = zip(gene_list, zipped_pvals)
    unpacked_gene_info = [(gene, pval, adjusted_pval) for (gene, (pval, adjusted_pval)) in zipped_gene_info]
    pval_df = pd.DataFrame(unpacked_gene_info, columns=['gene', 'pval', 'adjusted_pval'])
    print(len(pval_df), 'number of genes tested')

    pval_df.sort_values(by='pval', inplace=True)
    # put test name 
    pval_df.to_csv('allele_T-test_results_'+temp+'_MWU.txt', sep='\t', index=False)
    genes_to_graph_df = pval_df.loc[pval_df['adjusted_pval'] <= 0.05]  # set the pval cutoff for genes you want to make plots of
    genes_to_graph_list = genes_to_graph_df['gene'].tolist()


#RUN
if __name__ == "__main__":
    munged_files=argv[1:]
    for munged_file in munged_files:
        #print(munged_file)
        calc_pairs_and_compute_stats(munged_file)
        







