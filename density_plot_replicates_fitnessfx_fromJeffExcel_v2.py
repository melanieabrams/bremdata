import sys
import numpy as np

import seaborn as sns                                                
paper_rc = {'lines.linewidth': 10}                  
sns.set_context("paper", rc = paper_rc,font_scale=5)

import matplotlib.pyplot as plt

#usage: python3 density_plot_python.py rep_data.txt y.test_test_results

#Parameters

cutoff=False

cutoff_num = 10 #how many top genes to look at
specific_genes_to_include=['YGR098C','YMR168C','YHR023W','YLR397C']

def ParseFitness(rep_data):
    '''get the replicates' fitness data'''
    sp_dict={}
    sc_dict={}

    f = open(rep_data)
    next(f) #skip header row
    for line in f:
        #parse data
        line=line.strip()
        row_data=line.split("\t")
        #print(row_data)
        gene=row_data[-2]
        allele=row_data[-1]
        if allele=='sp':
            sp_dict[gene]=[x for x in row_data[1:-2] if x!='']
        elif allele=='sc':
            sc_dict[gene]=[x for x in row_data[1:-2] if x!='']

    f.close()

    return sp_dict,sc_dict


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


def PlotDensity(spar_ins, scer_ins,geneName):
    figName=geneName+"_py_density.pdf"

    fig=plt.figure(figsize=(20,15))

##    sns.distplot(spar_ins, hist = False, kde = True, kde_kws = {'linewidth':3},
##                 label= 'Insert in paradoxus')
##
##    sns.distplot(scer_ins, hist = False, kde = True, kde_kws = {'linewidth':3},
##                 label= 'Insert in cerevisiae')
    sns.kdeplot(spar_ins, bw='silverman', kernel='gau',label= 'Insert in paradoxus',color='#08A5CD')
    sns.kdeplot(scer_ins, bw='silverman', kernel='gau',label= 'Insert in cerevisiae',color='#F1B629')
    
    plt.title(geneName,fontsize=80,loc='center')
    plt.ylabel('Density of Replicates', fontsize=80)
    plt.xlabel('Fitness Score', fontsize=80)
    plt.rc('xtick',labelsize=80)
    plt.rc('ytick',labelsize=80)
    plt.savefig(figName, bbox_inches='tight', format='pdf',dpi=1000)

def PlotDensities(sp_dict,sc_dict,stats_dict=None,specific_genes=None):
    if stats_dict!=None:
        for gene in stats_dict:
            PlotDensity(sp_dict[gene], sc_dict[gene],gene)
    elif specific_genes!=None:
        for gene in specific_genes:
            if gene in sp_dict and gene in sc_dict:
                PlotDensity(sp_dict[gene], sc_dict[gene],gene)
            else:
                print(gene+" is not in both alleles")
    return 'done plotting'

rep_data=sys.argv[1]
sp_dict,sc_dict=ParseFitness(rep_data)

if cutoff!=False:
    stats_data=sys.argv[2]
    stats_dict=ParseResults(stats_data,cutoff_num)
    PlotDensities(sp_dict,sc_dict,stats_dict=stats_dict,specific_genes=specific_genes_to_include)

else:
    #print(sp_dict)
    #print(sc_dict['YHR023W'])
    #print(sp_dict['YHR023W'])
    PlotDensities(sp_dict,sc_dict,stats_dict=None,specific_genes=specific_genes_to_include)
    
