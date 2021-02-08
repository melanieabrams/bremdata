import sys
import numpy as np

import seaborn as sns                                                
paper_rc = {'lines.linewidth': 10}                  
sns.set_context("paper", rc = paper_rc,font_scale=5)

import matplotlib.pyplot as plt

#usage: python3 density_plot_python.py x.filtered_gene_inserts y.fxsize 

cutoff_num = 4 #how many top genes to look at

def ParseLogratios(insert_data,stats_dict):
    '''get the inserts' log2(39/28) for genes above cutoff'''
    sp_dict={}
    sc_dict={}
    for key in stats_dict:
        sp_dict[key]=[]
        sc_dict[key]=[]
    
    f = open(insert_data)
    header_line = True
    for line in f:
        #parse data
        row_data=line.split("\t")
        if header_line== True:
            gene_index = row_data.index('gene')
            allele_index = row_data.index('allele')
            if '39_28_log2' in row_data:  #this will be the header in collapsed pipeline
                logr_index=row_data.index('39_28_log2') 
            elif '39_28_log2_br_averaged_reads' in row_data: #this will be the header in uncollapsed
                logr_index=row_data.index('39_28_log2_br_averaged_reads')
            header_line = False
        else:
            gene=row_data[gene_index]
            if gene in stats_dict:
                allele=row_data[allele_index]
                logr=float(row_data[logr_index])
                #add to sc or sp dict
                if allele=='sp':
                    sp_dict[gene].append(logr)
                elif allele=='sc':
                    sc_dict[gene].append(logr)

    f.close()

    
    return sp_dict,sc_dict


def ParseResults(stats_data, cutoff_num):
    stats_dict={}
    with open (stats_data) as myfile:
        head = [next(myfile) for x in range(cutoff_num+1)]
        for line in head[1:]:
            row_data=line.split("\t")
            gene = row_data[0]
            fxsize=float(row_data[1].strip('\n'))
            stats_dict[gene]=fxsize

    return stats_dict


def PlotDensity(spar_ins, scer_ins,geneName):
    #figName=geneName+"_py_density.pdf"
    figName=geneName+"_py_density.png"

    fig=plt.figure(figsize=(20,15))

##    sns.distplot(spar_ins, hist = False, kde = True, kde_kws = {'linewidth':3},
##                 label= 'Insert in paradoxus')
##
##    sns.distplot(scer_ins, hist = False, kde = True, kde_kws = {'linewidth':3},
##                 label= 'Insert in cerevisiae')
    sns.kdeplot(spar_ins, bw='silverman', kernel='gau',label= 'Insert in paradoxus',color='#08A5CD')
    sns.kdeplot(scer_ins, bw='silverman', kernel='gau',label= 'Insert in cerevisiae',color='#F1B629')
    
    plt.title(geneName,fontsize=80,loc='center')
    plt.ylabel('Density of Inserts (kdeplot)', fontsize=80)
    plt.xlabel('log2(39/28)', fontsize=80)
    plt.xlim([-5,5])
    plt.rc('xtick',labelsize=80)
    plt.rc('ytick',labelsize=80)
    #plt.savefig(figName, bbox_inches='tight', format='pdf',dpi=1000)

    plt.savefig(figName, bbox_inches='tight', format='png')
    
def PlotDensities(stats_dict,sp_dict,sc_dict):
    for gene in stats_dict:
        PlotDensity(sp_dict[gene], sc_dict[gene],gene)
    return 'done plotting'

insert_data=sys.argv[1]
stats_data=sys.argv[2]

stats_dict=ParseResults(stats_data,cutoff_num)
sp_dict,sc_dict=ParseLogratios(insert_data,stats_dict)

PlotDensities(stats_dict,sp_dict,sc_dict)
