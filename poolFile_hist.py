import sys
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns                                                
paper_rc = {'lines.linewidth': 10}                  
sns.set_context("paper", rc = paper_rc,font_scale=5)

### USAGE ###
# python organize_and_filter_genes.py output_folder filtered_insert_ratio_file.poolcounts 
#adds rel_loc 
#splits tech reps and converts poolcount file formatted for this into .fastq_pooled_read_clean for downstream rh-seq pipeline

##PARAMETERS##
filename=sys.argv[1]
#filename='YBR061C_pf.txt'
maxCount=5000
figname='poolFile_bc_Hist.png'

##FUNCTIONS##

def parse_poolfile(filename, sep='\t',maxCount=5000):
    '''parses the poolcount file into a dataframe'''
    scer_ins=[]
    spar_ins=[]
    f=open(filename)
    next(f)
    for line in f:
        row_data=line.strip().split(sep)
        allele=row_data[4][:2]
        n=float(row_data[3])
        if n<=maxCount:
            if allele=='sc':
                scer_ins.append(n)
            elif allele=='sp':
                spar_ins.append(n)            
    f.close()
    #print(scer_ins,spar_ins)
    return scer_ins,spar_ins


def PlotHist(figName,scer_ins,spar_ins):
   # print(spar_ins)
    

    fig=plt.figure(figsize=(20,15))


    sns.kdeplot(spar_ins, bw='silverman', kernel='gau',label= 'Insert in paradoxus',color='#08A5CD')
    sns.kdeplot(scer_ins, bw='silverman', kernel='gau',label= 'Insert in cerevisiae',color='#F1B629')
    
    #plt.title(geneName,fontsize=80,loc='center')
    plt.ylabel('Density of BC Inserts', fontsize=80)
    plt.xlabel('Mean Insert Counts', fontsize=80)
    plt.rc('xtick',labelsize=80)
    plt.rc('ytick',labelsize=80)
    #plt.xscale('log')

    #plt.vlines(0,0,1,colors='grey',linestyles='dashed',linewidth=5)

    #plt.xlim(-5,5)
    plt.savefig(figName, bbox_inches='tight', format='png')
    print("done plotting histogram: "+figName)

        

##START
if __name__ == '__main__':
    print('...parsing poolcount...')
    scer_ins,spar_ins=parse_poolfile(filename,maxCount=10000000000000)
    print('...plotting...')
    PlotHist(figname,scer_ins,spar_ins)
