import sys
import numpy as np

from itertools import chain

import seaborn as sns                                                
paper_rc = {'lines.linewidth': 10}                  
sns.set_context("paper", rc = paper_rc,font_scale=5)

import matplotlib.pyplot as plt

#usage: python3 moreplots.py

#Parameters


fgi='2.0_10_coeffvarANDreadcutoff_uncollapsedMBA_unpaired_ratiocvfilt_10.0_2.0.filtered_gene_inserts'
figname='2pt0_10_10_2_fgi_38C_hist_atleast9br'
temp='38'
analysis='RBseq_to_RHseq'

degree_sign='\u00b0'

#Functions



def ParseFGI(filtered_gene_inserts,rep,analysis=None):
    '''get the inserts' log2(39/28)'''
    sp_dict={}
    sc_dict={}

    print(filtered_gene_inserts,analysis)

    YAL002W_ins=0
    print("YAL002W sc insert values: ")
    f = open(filtered_gene_inserts)
    next(f) #skip header row
    for line in f:
        #parse data
        line=line.strip()
        row_data=line.split("\t")

        if analysis=='RBseq_to_RHseq':
            gene=row_data[-3]
            allele=row_data[-2]
            ins_ID=row_data[-18]
            logratio_avg_39to28=row_data[-18+rep]
            br12=row_data[-6]
            reads39C=row_data[-41+rep]
            brs_present=len([x for x in row_data[-41:-29] if x!=''])
            if logratio_avg_39to28!='' and brs_present>=9:
            #if logratio_avg_39to28!='' and br12!='':
            #if logratio_avg_39to28!='':
                logratio_avg_39to28=float(logratio_avg_39to28)
                if float(reads39C)>=1:
                    if gene=="YAL002W" and allele=="sc":
                        print(ins_ID,logratio_avg_39to28)
                        YAL002W_ins+=1
                    if allele=='sp':
                        sp_dict[ins_ID]=logratio_avg_39to28
                    elif allele=='sc':
                        sc_dict[ins_ID]=logratio_avg_39to28
                    


##
    print('YAL002W ins: ',YAL002W_ins)
    print('mean of sc inserts: ',sum(sc_dict.values())/len(sc_dict))
    
    print('num of sc inserts', len(sc_dict))

    print('mean of sp inserts: ',sum(sp_dict.values())/len(sp_dict))
    
    print('num of sp inserts', len(sp_dict))

    print('num_ins_total',(len(sc_dict)+len(sp_dict)))
    
##    print(sum(sc_dict.values()))

    f.close()

    return sp_dict,sc_dict



def PlotHist(figName,sp_dict,sc_dict, rep,onlyGOI=False):

    fig=plt.figure(figsize=(20,15))

    spar_ins=[]
    for key in sp_dict:
        if onlyGOI==False:
            spar_ins.append(sp_dict[key])
        elif key in genes_to_bar:
            spar_ins.append(sp_dict[key])
    scer_ins=[]
    for key in sc_dict:
        if onlyGOI==False:
            scer_ins.append(sc_dict[key])
        elif key in genes_to_bar:
            scer_ins.append(sc_dict[key])

    sns.kdeplot(spar_ins, bw='silverman', kernel='gau',label= 'Insert in paradoxus',color='#08A5CD')
    sns.kdeplot(scer_ins, bw='silverman', kernel='gau',label= 'Insert in cerevisiae',color='#F1B629')
    
    #plt.title(geneName,fontsize=80,loc='center')
    plt.ylabel('Density of Inserts', fontsize=80)
    plt.xlabel('logs('+temp+'/28'+degree_sign+'C) Averaged Reads)', fontsize=80)
    plt.rc('xtick',labelsize=80)
    plt.rc('ytick',labelsize=80)

    plt.vlines(0,0,1,colors='grey',linestyles='dashed',linewidth=5)

    plt.savefig(figName+'rep_'+str(i)+'.png', bbox_inches='tight', format='png')
    print("done plotting histogram: "+figName+'rep_'+str(rep)+'.png')




#RUN


for i in range(1,13):
    #parse data
    sp_dict,sc_dict=ParseFGI(fgi,i,analysis=analysis)
    #plot hist
    PlotHist(figname,sp_dict,sc_dict,i)

