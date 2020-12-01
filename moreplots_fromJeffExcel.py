import sys
import numpy as np

import seaborn as sns                                                
paper_rc = {'lines.linewidth': 10}                  
sns.set_context("paper", rc = paper_rc,font_scale=5)

import matplotlib.pyplot as plt

#usage: python3 moreplots.py

data_36='geneFit_36C_dropNA_allele_sorted.csv'
data_37='geneFit_37C_dropNA_allele_sorted.csv'
figName='36_and_37_MeansStdevScatter.pdf'

#Functions

def ParseFitness(rep_data):
    '''get the inserts' log2(39/28) for genes above cutoff'''
    sp_dict={}
    sc_dict={}

    f = open(rep_data)
    next(f) #skip header row
    for line in f:
        #parse data
        line=line.strip()
        row_data=line.split(",")
        gene=row_data[-2]
        allele=row_data[-1]
        rep_vals=[]
        for i in row_data[1:-2]:
            rep_vals.append(float(i))
        if allele=='sp':
            sp_dict[gene]=rep_vals
        elif allele=='sc':
            sc_dict[gene]=rep_vals

    f.close()

    return sp_dict,sc_dict

def calcTStats(sp_dict, sc_dict):
    sp_T_dict={}
    sc_T_dict={}

    #get mean number and size of observations across everything
    samplesize=0
    sumvals=0
    for gene in sp_dict:
        sumvals+=np.sum(sp_dict[gene])
        samplesize+=len(sp_dict[gene])
    for gene in sc_dict:
        sumvals+=np.sum(sc_dict[gene])
        samplesize+=len(sc_dict[gene])
    sampleavg=float(sumvals)/float(samplesize)
        

    #calc t stat
    
    for gene in sp_dict:
        reps=sp_dict[gene]
        T=np.mean(reps)/np.std(reps)  #NOTE THAT THIS IS NOT WHAT JEFF HAS NOW - IS THIS THE T-stat
        #T=(np.mean(reps)-sampleavg)/(np.std(reps)/np.sqrt(samplesize))#  THIS IS ALSO NOT RIGHT, AS I TRIED to GUESS WHAT HE WOULD HAVE DONE
        sp_T_dict[gene]=T

    for gene in sc_dict:
        reps=sc_dict[gene]
        T=np.mean(reps)/np.std(reps)
        sc_T_dict[gene]=T

    #print(sp_T_dict['YLR397C'])

    return sp_T_dict,sc_T_dict

def calcRatios(sp_T_dict,sc_T_dict):
    '''for all genes with both sc and sp, calc Tstat ratios'''
    return Tratio_dict


def PlotScatter37_36(figName,sp_T_dict_36,sc_T_dict_36,sp_T_dict_37,sc_T_dict_37):

    fig=plt.figure(figsize=(20,15))

    x_sp=[]
    y_sp=[]
    for key in sp_T_dict_36:
        if key in sp_T_dict_37:
            x_sp.append(sp_T_dict_36[key])
            y_sp.append(sp_T_dict_37[key])

    x_sc=[]
    y_sc=[]
    for key in sc_T_dict_36:
        if key in sc_T_dict_37:
            x_sc.append(sc_T_dict_36[key])
            y_sc.append(sc_T_dict_37[key])

    line1=plt.scatter(x_sp, y_sp, color='#F1B629')
    line2=plt.scatter(x_sc, y_sc, color='#08A5CD')

    plt.hlines(0,-4,4,colors='grey',linestyles='dashed', linewidth=1)
    plt.vlines(0,-4,4,colors='grey',linestyles='dashed',linewidth=1)

    plt.legend((line1, line2), ('sp', 'sc'))
    
    plt.title(figName.split('.')[0],fontsize=80,loc='center')
    plt.ylabel('37 means/sdtevs', fontsize=80)
    plt.xlabel('36 means/sdtevs', fontsize=80)
    plt.rc('xtick',labelsize=80)
    plt.rc('ytick',labelsize=80)
    plt.savefig(figName, bbox_inches='tight', format='pdf',dpi=1000)




###RUN

sp_dict_37,sc_dict_37=ParseFitness(data_37)
sp_T_dict_37,sc_T_dict_37=calcTStats(sp_dict_37,sc_dict_37)

sp_dict_36,sc_dict_36=ParseFitness(data_36)
sp_T_dict_36,sc_T_dict_36=calcTStats(sp_dict_36,sc_dict_36)
               
PlotScatter37_36(figName,sp_T_dict_36,sc_T_dict_36,sp_T_dict_37,sc_T_dict_37)
