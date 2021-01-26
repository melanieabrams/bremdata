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
#filename='YBR061C_poolCount_cvfilt3pt3.txt'
maxCount=20000000
figname='poolCount_Hist_maxCt'+str(maxCount)+'.png'

##FUNCTIONS##

def parse_poolcount(filename, sep='\t',maxCount=1000):
    '''parses the poolcount, looks for the maximum value'''
    scer_ins=[]
    spar_ins=[]
    f=open(filename)
    next(f)
    for line in f:
        row_data=line.strip().split(sep)
        allele=row_data[1][:2]
        brs=list(map(float, row_data[8:]))
        #print(brs)
        n=np.max(brs)
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
    plt.ylabel('Density of Inserts', fontsize=80)
    plt.xlabel('PoolCount Max Count for Ins', fontsize=80)
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
    scer_ins,spar_ins=parse_poolcount(filename,maxCount=maxCount)
    print('...plotting...')
    PlotHist(figname,scer_ins,spar_ins)
