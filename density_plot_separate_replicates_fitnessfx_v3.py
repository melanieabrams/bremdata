import sys
import numpy as np
import copy

import pandas as pd

import seaborn as sns                                                
paper_rc = {'lines.linewidth': 10}                  
#sns.set_context("paper", rc = paper_rc,font_scale=5)

import matplotlib.pyplot as plt

#usage: python3 density_plot_python.py rep_data.txt y.test_test_results

#Parameters

rep_data='37C_2.0_1.1_coeffvarANDreadcutoff_uncollapsedMBA_unpaired_ratiocvfilt_10.0_3.0.filtered_gene_inserts'
stats_data='37C_2.0_1.1_coeffvarANDreadcutoff_uncollapsedMBA_unpaired_ratiocvfilt_10.0_3.0.mwu_test_results'

cutoff=False

cutoff_num = 10 #how many top genes to look at
specific_genes_to_include=['YLR397C']


def ParseFitness(rep_data):
    '''get the replicates' fitness data'''
    sp_dict={}
    sc_dict={}

    df=pd.read_csv(rep_data,sep='\t')

    print(df.head())
    print(df.columns)

    fitness_ratio_rep_columns=[col for col in df.columns if '_n_av' in col and '39_28_log2' in col] # get columns with reps' data
    #initialize each rep in dict for each allele:
    empty_rep_dict={}
    for col in fitness_ratio_rep_columns:
        empty_rep_dict[col]=[]
    all_genes=list(set(df['gene'].to_list()))
    for gene in all_genes:
        sc_dict[gene]=copy.deepcopy(empty_rep_dict)
        sp_dict[gene]=copy.deepcopy(empty_rep_dict)

    #add all observations

    for index, row in df.iterrows():
        allele=row['allele']
        gene=row['gene']
        if allele=='sc':
            for rep_col in fitness_ratio_rep_columns:
                sc_dict[gene][rep_col].append(row[rep_col])
      

        elif allele=='sp':
            for rep_col in fitness_ratio_rep_columns: 
                 sp_dict[gene][rep_col].append(row[rep_col])
    

        

    return sp_dict,sc_dict,fitness_ratio_rep_columns


def ParseResults(stats_data, cutoff_num):
    stats_dict={}
    with open (stats_data) as myfile:
        #get top # of genes
        head = [next(myfile) for x in range(cutoff_num+1)]
        for line in head[1:]:
            row_data=line.split("\t")
            gene = row_data[0]
            adj_pval=float(row_data[2].strip('\n'))
            stats_dict[gene]=adj_pval
        #get specific genes
        for line in myfile:
            row_data=line.split("\t")
            gene = row_data[0]
            adj_pval=float(row_data[2].strip('\n'))
            if gene in specific_genes_to_include:
                stats_dict[gene]=adj_pval
    return stats_dict


def PlotDensity(spar_ins, scer_ins,rep_columns,geneName):
    figName=geneName+"_py_density_separate_replicates.pdf"

    fig, axs = plt.subplots(3,4, figsize=(20, 15))
    fig.subplots_adjust(hspace = .5, wspace=.5)
    fig.suptitle(geneName)
    plt.rc('xtick',labelsize=10)
    plt.rc('ytick',labelsize=10)


    for ax, rep_col in zip(axs.ravel(), rep_columns):
            

        sns.kdeplot(spar_ins[rep_col], ax=ax,bw_method='silverman',label= 'Insert in paradoxus',color='#08A5CD') #only gaussian kernels supported, but that's what I'd been using anyway
        sns.kdeplot(scer_ins[rep_col], ax=ax,bw_method='silverman',label= 'Insert in cerevisiae',color='#F1B629') #also bw_method replacing deprecated bw
        
        ax.set_title('Replicate '+rep_col[10])
        ax.set_ylabel('Density of Observations', fontsize=10)
        ax.set_xlabel('Fitness Score', fontsize=10)
            
    #plt.show()
    plt.savefig(figName, bbox_inches='tight', format='pdf',dpi=1000)


def PlotDensities(sp_dict,sc_dict,rep_columns,stats_dict=None,specific_genes=None):
    if stats_dict!=None:
        for gene in stats_dict:
            PlotDensity(sp_dict[gene], sc_dict[gene],rep_columns,gene)
    elif specific_genes!=None:
        for gene in specific_genes:
            if gene in sp_dict and gene in sc_dict:
                PlotDensity(sp_dict[gene], sc_dict[gene],rep_columns,gene)
            else:
                print(gene+" is not in both alleles")
    return 'done plotting'


##RUN###
sp_dict,sc_dict,rep_columns=ParseFitness(rep_data)
print(sp_dict)

if cutoff!=False:
    stats_dict=ParseResults(stats_data,cutoff_num)
    PlotDensities(sp_dict,sc_dict,rep_columns,stats_dict=stats_dict,specific_genes=specific_genes_to_include)

else:
    #print(sp_dict)
    #print(sc_dict['YHR023W'])
    #print(sp_dict['YHR023W'])
    PlotDensities(sp_dict,sc_dict,rep_columns,stats_dict=None,specific_genes=specific_genes_to_include)
    
