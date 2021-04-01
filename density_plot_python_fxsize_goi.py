import sys
import numpy as np

import seaborn as sns                                                
paper_rc = {'lines.linewidth': 10}                  
sns.set_context("paper", rc = paper_rc,font_scale=5)

import matplotlib.pyplot as plt

#usage: python3 density_plot_python.py x.filtered_gene_inserts y.fxsize 



temperature='37'
#goi: 2_1.1_1.6_3.0, pcutoff 0.05, abs(sp_mean_obs)<=0.75, sc_mean_obs<-0.5
goi=['YGR198W','YLR397C','YBR136W','YGR140W','YNL049C','YKL017C','YGL082W',
     'YMR125W','YDR508C','YMR207C','YDR375C','YDR180W','YKL197C','YDR318W',
     'YMR094W','YOR326W','YBR081C','YPR049C','YIL152W','YER151C','YJR107W',
     'YAL026C','YDR456W','YLR141W','YPL268W','YDR235W']

cutoff_num=5 #use to only look at top # of goi
#cutoff_num=len(goi) #uncomment to look at all goi

insert_data='37C_2.0_1.1_coeffvarANDreadcutoff_uncollapsedMBA_unpaired_ratiocvfilt_1.6_3.0.filtered_gene_inserts'
stats_data='37C_2.0_1.1_coeffvarANDreadcutoff_uncollapsedMBA_unpaired_ratiocvfilt_1.6_3.0.all_fxsizes'

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


def ParseResults(stats_data,goi,cutoff_num):
    stats_dict={}
    with open (stats_data) as myfile:
        for line in myfile:
            row_data=line.split("\t")
            gene = row_data[0]
            if gene in goi:
                fxsize=float(row_data[1].strip('\n'))
                stats_dict[gene]=fxsize
                if len(stats_dict)==cutoff_num:
                    return stats_dict
    return stats_dict


def PlotDensity(spar_ins, scer_ins,geneName):
    #figName=geneName+"_py_density.pdf"
    figName=temperature+'_'+geneName+"_py_density.png"

    fig=plt.figure(figsize=(20,15))

    sns.kdeplot(spar_ins,label= 'Insert in paradoxus',color='#08A5CD')
    sns.kdeplot(scer_ins,label= 'Insert in cerevisiae',color='#F1B629')
    
    plt.title(geneName,fontsize=80,loc='center')
    plt.ylabel('Density of Inserts', fontsize=10)
    plt.xlabel('log2('+temperature+'/28)', fontsize=10)
    plt.xlim([-5,5])
    plt.rc('xtick',labelsize=80)
    plt.rc('ytick',labelsize=80)
    #plt.savefig(figName, bbox_inches='tight', format='pdf',dpi=1000)

    plt.savefig(figName, bbox_inches='tight', format='png')
    
def PlotDensities(stats_dict,sp_dict,sc_dict):
    for gene in stats_dict:
        print(gene)
        PlotDensity(sp_dict[gene], sc_dict[gene],gene)
    return 'done plotting'

###RUN###

stats_dict=ParseResults(stats_data,goi,cutoff_num)
sp_dict,sc_dict=ParseLogratios(insert_data,stats_dict)

PlotDensities(stats_dict,sp_dict,sc_dict)
