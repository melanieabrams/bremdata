import sys
import numpy as np

import seaborn as sns                                                
paper_rc = {'lines.linewidth': 10}                  
sns.set_context("paper", rc = paper_rc,font_scale=5)

import matplotlib.pyplot as plt

#usage: python3 density_plot_python.py x.filtered_gene_inserts y.fxsize 



temperature='39'
#goi: 2_1.1_1.6_3.0, pcutoff 0.05, abs(sp_mean_obs)<=0.75, sc_mean_obs<-0.5
goi=['YGR198W', 'YMR207C', 'YGL082W', 'YNL049C', 'YDL035C', 'YDR508C', 'YBR136W',
 'YML099C', 'YPL254W', 'YIL152W', 'YKL017C', 'YGR140W', 'YJR127C', 'YDR375C',
 'YOR091W', 'YLR397C', 'YNL132W', 'YMR078C', 'YLR422W', 'YMR125W', 'YOR371C',
 'YMR094W', 'YMR167W', 'YDR103W', 'YDR318W', 'YAL026C', 'YDR180W', 'YOR092W',
 'YDR235W', 'YER151C', 'YMR275C', 'YKL114C', 'YOL081W', 'YPR049C', 'YGL095C',
 'YDR456W', 'YKL197C', 'YIL068C', 'YOR326W', 'YNR045W', 'YJR107W', 'YPL268W',
 'YJL062W', 'YCR042C'] #hits 37_2.0_1.1_..._10.0_3.0 with sc defect and effect size >0.5

cutoff_num=44 #use to only look at top # of goi
#cutoff_num=len(goi) #uncomment to look at all goi

insert_data='Carly_Reanalysis_1.5_1.1_coeffvarANDreadcutoff_uncollapsedMBA_unpaired_ratiocvfilt_20.0_2.0.filtered_gene_inserts'
stats_data='Carly_Reanalysis_1.5_1.1_coeffvarANDreadcutoff_uncollapsedMBA_unpaired_ratiocvfilt_20.0_2.0.all_fxsizes'

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
