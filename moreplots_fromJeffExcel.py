import sys
import numpy as np

import seaborn as sns                                                
paper_rc = {'lines.linewidth': 10}                  
sns.set_context("paper", rc = paper_rc,font_scale=5)

import matplotlib.pyplot as plt

#usage: python3 moreplots.py

#Parameters

data_36='geneFit_36C_dropNA_allele_sorted.csv'
gen_36=[12.84584039,12.96574443,12.82881066,12.47694361,12.69612127,12.57667477,12.72712914,12.42038924,12.4493805,12.52735329,13.08141686,13.31136996]

data_37='geneFit_37C_dropNA_allele_sorted.csv'
gen_37=[10.71381697,10.6059797,10.70591248,10.74969728,10.72307974,10.87017294,10.95210186,10.61683005,9.636340693,10.18646952,10.72492607,10.55877579]

scatterFigName='36_and_37_MeansStdevScatter.pdf'

degree_sign='\u00b0'

barFigName='36_and_37_fitnessBars.png'
genes_to_bar=['YGR098C','YMR168C','YHR023W','YLR397C','YKR054C','YDR180W','YGR198W','YHR166C','YCR042C','YPL174C','YEL030W','YOR151C','YNL027W','YBR136W','YGL205W']


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


def normToGen(rep_dict,gen):
    norm={}

    for gene in rep_dict:
        norm[gene]=[]
        for i in range(len(rep_dict[gene])):
            norm[gene].append(float(rep_dict[gene][i])/gen[i])

    return norm
            



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

def PlotBars(figName,dict_36,dict_37,norm=False):
    '''bar plots of genes of interest' fitness for the allele passed in the dict''' #note: modified from stationary/active plot, hence A and S naming

    w=0.35 #bar width

    fontsize=16

    colors=[] #replace colors with B/W
    for i in range(len(genes_to_bar)):
        colors.append((0,0,0))

    #colors for replicates' dots
    A_colors=[] 
    for i in range(len(genes_to_bar)):
        A_colors.append((0.5, 0.75, 1, 1))

    S_colors=[] 
    for i in range(len(genes_to_bar)):
        S_colors.append((0.75, 0.75, 0.75, 1))
    
        
    x=range(len(genes_to_bar)) #x: num genes of interest
    yS=[] #y: mean of bioreps for the strains, A for 36, S for 37
    yA=[]
    Alabels= [ ] #label: gene names in order
    Slabels= [ ]


    for gene in genes_to_bar:
        yA.append(np.array(dict_36[gene]))
        Alabels.append(gene)
    for gene in genes_to_bar:
        yS.append(np.array(dict_37[gene]))
        Slabels.append(gene)

          
    fig, ax = plt.subplots(figsize=(18,6))
    ind=np.arange(len(yA))

    rectsA=ax.bar(ind-w/2,
       height=[np.mean(yi) for yi in yA],
       width=w,    # bar width
       label='36'+degree_sign, color=(0.5, 0.75, 1, 0.75),  
       edgecolor=colors,zorder=-1
       )

    rectsS=ax.bar(ind+w/2,
       height=[np.mean(yi) for yi in yS],
       width=w, label='37'+degree_sign, color=(0.75, 0.75, 0.75, 0.75), 
       edgecolor=colors,zorder=-1
       )

 

    ax.set_xticks(ind)
    ax.set_xticklabels(genes_to_bar,fontsize=8)
    ax.legend(fontsize=8)

    plt.yticks(fontsize=8)
    plt.ylim(-10,10)
    ax.set_yticks(range(-10,10))
    if norm==True:
        plt.ylim(-1,1)
        ax.set_yticks(range(-1,1))

    #ax.set_yticks(range(-5,5))


    #scatter individual strains yA
    for i in range(int(len(x))):
        nonzero_yi=np.where(yA[i]==0, 1, yA[i]) 
        yA_scatter=ax.scatter(x[i] + np.full_like(nonzero_yi,0.5) * w - w / 2 - w/2, nonzero_yi, color=A_colors[i],zorder=1)

    #scatter individual strains yS
    for i in range(int(len(x))):
        nonzero_yi=np.where(yS[i]==0, 1, yS[i]) 
        yS_scatter=ax.scatter(x[i] + np.full_like(nonzero_yi,0.5) * w - w / 2 + w/2, nonzero_yi, color=S_colors[i],zorder=1)

   
    x_for_line_at_1 = [ax.patches[0].get_x(), ax.patches[-1].get_x() + ax.patches[-1].get_width()] #https://stackoverflow.com/questions/38017465/how-to-add-a-line-on-top-of-a-bar-chart
    ax.plot(x_for_line_at_1, [0,0], 'r--', c='grey', linewidth=1)

    plt.ylabel('Raw Fitness Score (Temp vs. 28'+degree_sign+'C)', fontsize=fontsize)
    plt.xlabel('GOI', fontsize=fontsize)

##    plt.yscale('log')
##    plt.ylim(10**2,5*10**7)

    print('done plotting bar plot: '+figName)
    plt.tight_layout()

    plt.savefig(figName)
    
    return      



###RUN

sp_dict_37,sc_dict_37=ParseFitness(data_37)
sp_T_dict_37,sc_T_dict_37=calcTStats(sp_dict_37,sc_dict_37)
sp_norm_dict_37,sc_norm_dict_37=normToGen(sp_dict_37,gen_37),normToGen(sc_dict_37,gen_37)

sp_dict_36,sc_dict_36=ParseFitness(data_36)
sp_T_dict_36,sc_T_dict_36=calcTStats(sp_dict_36,sc_dict_36)
sp_norm_dict_36,sc_norm_dict_36=normToGen(sp_dict_36,gen_36),normToGen(sc_dict_36,gen_36)

####scatterplot of mean/stdev at 36 vs 37               
#PlotScatter37_36(scatterFigName,sp_T_dict_36,sc_T_dict_36,sp_T_dict_37,sc_T_dict_37)

####barplots: unnormalized
##PlotBars("sp_"+barFigName,sp_dict_36,sp_dict_37)
##PlotBars("sc_"+barFigName,sc_dict_36,sc_dict_37)

####barplots: fitness/generations
PlotBars("sp_normToGen_"+barFigName,sp_norm_dict_36,sp_norm_dict_37,norm=True)
PlotBars("sc_normToGen"+barFigName,sc_norm_dict_36,sc_norm_dict_37,norm=True)
