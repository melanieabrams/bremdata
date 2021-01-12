import sys
import numpy as np

from itertools import chain

import seaborn as sns                                                
paper_rc = {'lines.linewidth': 10}                  
sns.set_context("paper", rc = paper_rc,font_scale=5)

import matplotlib.pyplot as plt

#usage: python3 moreplots.py

#Parameters

data_36='geneFit_36C_dropNA.txt'
gen_36=[12.84584039,12.96574443,12.82881066,12.47694361,12.69612127,12.57667477,12.72712914,12.42038924,12.4493805,12.52735329,13.08141686,13.31136996]

data_37='geneFit_37C_dropNA.txt'
gen_37=[10.71381697,10.6059797,10.70591248,10.74969728,10.72307974,10.87017294,10.95210186,10.61683005,9.636340693,10.18646952,10.72492607,10.55877579]

data_38='geneFit_38C_dropNA.txt'
gen_38=[5.123086751,5.10433666,5.121015401,5.203592714,5.037821465,5.022367813,4.83541884,4.882643049,4.832890014,4.794415866,4.860466259,4.87036472]



scatterFigName='36_and_37_MeansStdevScatter.png'
scatterFigName_for38='37_and_38_MeanStdevScatter.png'

degree_sign='\u00b0'

barFigName='36_and_37_fitnessBars.png'
barFigName3='36_37_and_38_fitnessBars.png'

genes_to_bar=['YGR098C','YMR168C','YHR023W','YLR397C','YKR054C','YDR180W','YGR198W','YHR166C','YCR042C','YPL174C','YOR151C','YBR136W','YGL205W']

#genes not in dataset (removed from plots): 'YEL030W','YNL027W'


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
        row_data=line.split("\t")
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

def calcMean(sp_dict, sc_dict):
    sp_M_dict={}
    sc_M_dict={}
    
    for gene in sp_dict:
        reps=sp_dict[gene]
        M=np.mean(reps)
        sp_M_dict[gene]=M

    for gene in sc_dict:
        reps=sc_dict[gene]
        M=np.mean(reps)
        sc_M_dict[gene]=M


    return sp_M_dict,sc_M_dict


def normToGen(rep_dict,gen):
    norm={}

    for gene in rep_dict:
        norm[gene]=[]
        for i in range(len(rep_dict[gene])):
            norm[gene].append(float(rep_dict[gene][i])/gen[i])

    return norm

def normToMean(rep_dict):
    norm={}

    rep_fit_vals=[]
    for gene in rep_dict:
        rep_fit_vals+=rep_dict[gene]

    meanFitness=np.mean(rep_fit_vals)

    for gene in rep_dict:
        norm[gene]=[]
        for i in range(len(rep_dict[gene])):
            norm[gene].append(float(rep_dict[gene][i])/meanFitness)

    return norm
            


def PlotScatter37_36(figName,sp_T_dict_36,sc_T_dict_36,sp_T_dict_37,sc_T_dict_37,onlyGOI=False,justMean=False):

    fig=plt.figure(figsize=(20,15))

    x_sp=[]
    y_sp=[]
    for key in sp_T_dict_36:
        if key in sp_T_dict_37:
            if onlyGOI==False:
                x_sp.append(sp_T_dict_36[key])
                y_sp.append(sp_T_dict_37[key])
            elif key in genes_to_bar:
                x_sp.append(sp_T_dict_36[key])
                y_sp.append(sp_T_dict_37[key])

    x_sc=[]
    y_sc=[]
    for key in sc_T_dict_36:
        if key in sc_T_dict_37:
            if onlyGOI==False:
                x_sc.append(sc_T_dict_36[key])
                y_sc.append(sc_T_dict_37[key])
            elif key in genes_to_bar:
                x_sc.append(sc_T_dict_36[key])
                y_sc.append(sc_T_dict_37[key])
            
    line1=plt.scatter(x_sp, y_sp, color='#F1B629')
    line2=plt.scatter(x_sc, y_sc, color='#08A5CD')

    plt.hlines(0,-4,4,colors='grey',linestyles='dashed', linewidth=1)
    plt.vlines(0,-4,4,colors='grey',linestyles='dashed',linewidth=1)

    plt.legend((line1, line2), ('sp', 'sc'))
    
    
    if justMean==True:
        plt.ylabel('37 Mean Fitness', fontsize=80)
        plt.xlabel('36 Mean Fitness', fontsize=80)
        #plt.title(figName.split('Stdev')[0],fontsize=80,loc='center')

    else:
        #plt.title(figName.split('.')[0],fontsize=80,loc='center')
        plt.ylabel('37 means/sdtevs', fontsize=80)
        plt.xlabel('36 means/sdtevs', fontsize=80)
        
    plt.rc('xtick',labelsize=80)
    plt.rc('ytick',labelsize=80)
    plt.savefig(figName, bbox_inches='tight', format='png')
    print("done plotting scatterplot: "+figName)


def PlotScatter37_38(figName,sp_T_dict_37,sc_T_dict_37,sp_T_dict_38,sc_T_dict_38,onlyGOI=False,justMean=False):

    fig=plt.figure(figsize=(20,15))

    x_sp=[]
    y_sp=[]
    for key in sp_T_dict_37:
        if key in sp_T_dict_38:
            if onlyGOI==False:
                x_sp.append(sp_T_dict_37[key])
                y_sp.append(sp_T_dict_38[key])
            elif key in genes_to_bar:
                x_sp.append(sp_T_dict_37[key])
                y_sp.append(sp_T_dict_38[key])

    x_sc=[]
    y_sc=[]
    for key in sc_T_dict_37:
        if key in sc_T_dict_38:
            if onlyGOI==False:
                x_sc.append(sc_T_dict_37[key])
                y_sc.append(sc_T_dict_38[key])
            elif key in genes_to_bar:
                x_sc.append(sc_T_dict_37[key])
                y_sc.append(sc_T_dict_38[key])
            
    line1=plt.scatter(x_sp, y_sp, color='#08A5CD')
    line2=plt.scatter(x_sc, y_sc, color='#FF8000')

    plt.hlines(0,-4,4,colors='grey',linestyles='dashed', linewidth=1)
    plt.vlines(0,-4,4,colors='grey',linestyles='dashed',linewidth=1)

    plt.legend((line1, line2), ('sp', 'sc'), prop={'size': 80})
    
    #plt.title(figName.split('.')[0],fontsize=80,loc='center')
    if justMean==True:
        plt.ylabel('38 Mean Fitness', fontsize=80)
        plt.xlabel('37 Mean Fitness', fontsize=80)
    else:
        plt.ylabel('38 means/sdtevs', fontsize=80)
        plt.xlabel('37 means/sdtevs', fontsize=80)
    plt.rc('xtick',labelsize=80)
    plt.rc('ytick',labelsize=80)
    plt.savefig(figName, bbox_inches='tight', format='png')
    print("done plotting scatterplot: "+figName)



def Plot2Bars(figName,dict_36,dict_37,norm=False):
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
       width=w, label='37'+degree_sign+'C', color=(0.75, 0.75, 0.75, 0.75), 
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

    if norm==True:
        plt.ylabel('Generation-Normalized Fitness Score (Temp vs. 28'+degree_sign+'C)', fontsize=fontsize)
    else:
        plt.ylabel('Raw Fitness Score (Temp vs. 28'+degree_sign+'C)', fontsize=fontsize)
    plt.xlabel('GOI', fontsize=fontsize)

##    plt.yscale('log')
##    plt.ylim(10**2,5*10**7)

    print('done plotting bar plot: '+figName)
    plt.tight_layout()

    plt.savefig(figName)
    
    return


def Plot3Bars(figName,dict_36,dict_37,dict_38,normToGen=False,normToMean=False):
    '''bar plots of genes of interest' fitness for the allele passed in the dict''' #note: modified from stationary/active plot, hence A and S naming
        

    w=0.2 #bar width

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

    H_colors=[]
    for i in range(len(genes_to_bar)):
        H_colors.append((1,0.5,0,1))
    
        
    x=range(len(genes_to_bar)) #x: num genes of interest
    yS=[] #y: mean of bioreps for the strains, A for 36, S for 37
    yA=[]
    yH=[]
    Alabels= [ ] #label: gene names in order
    Slabels= [ ]
    Hlabels= [ ]


    for gene in genes_to_bar:
        yA.append(np.array(dict_36[gene]))
        Alabels.append(gene)
    for gene in genes_to_bar:
        yS.append(np.array(dict_37[gene]))
        Slabels.append(gene)
    for gene in genes_to_bar:
        yH.append(np.array(dict_38[gene]))
        Hlabels.append(gene)

          
    fig, ax = plt.subplots(figsize=(18,6))
    ind=np.arange(len(yA))

    rectsA=ax.bar(ind-w,
       height=[np.mean(yi) for yi in yA],
       width=w,    # bar width
       label='36'+degree_sign, color=(0.5, 0.75, 1, 0.75),  
       edgecolor=colors,zorder=-1
       )

    rectsS=ax.bar(ind+0,
       height=[np.mean(yi) for yi in yS],
       width=w, label='37'+degree_sign+'C', color=(0.75, 0.75, 0.75, 0.75), 
       edgecolor=colors,zorder=-1
       )

    rectsH=ax.bar(ind+w,
       height=[np.mean(yi) for yi in yH],
       width=w, label='38'+degree_sign+'C', color=(1.0, 0.5, 0.0, 0.75), 
       edgecolor=colors,zorder=-1
       )
 

    ax.set_xticks(ind)
    ax.set_xticklabels(genes_to_bar,fontsize=8)
    ax.legend(fontsize=8)

    
 
    #scatter individual strains yA
    for i in range(int(len(x))):
        nonzero_yi=np.where(yA[i]==0, 1, yA[i]) 
        yA_scatter=ax.scatter(x[i] + np.full_like(nonzero_yi,0.5) * w - w / 2 - w, nonzero_yi, color=A_colors[i],zorder=1)

    #scatter individual strains yS
    for i in range(int(len(x))):
        nonzero_yi=np.where(yS[i]==0, 1, yS[i]) 
        yS_scatter=ax.scatter(x[i] + np.full_like(nonzero_yi,0.5) * w - w / 2 +0 , nonzero_yi, color=S_colors[i],zorder=1)

    #scatter individual strains yH
    for i in range(int(len(x))):
        nonzero_yi=np.where(yH[i]==0, 1, yH[i]) 
        yH_scatter=ax.scatter(x[i] + np.full_like(nonzero_yi,0.5) * w - w / 2 + w , nonzero_yi, color=H_colors[i],zorder=1)
   
    x_for_line_at_1 = [ax.patches[0].get_x(), ax.patches[-1].get_x() + ax.patches[-1].get_width()] #https://stackoverflow.com/questions/38017465/how-to-add-a-line-on-top-of-a-bar-chart
    ax.plot(x_for_line_at_1, [0,0], 'r--', c='grey', linewidth=1)

    #get and set  plottable ylims
    maxY=(max(chain(max(i) for i in (yA+yH+yS))))
    minY=(min(chain(min(i) for i in (yA+yH+yS))))
    bound=max(abs(maxY),abs(minY))
    bound=bound+round(bound/10,1)
    plt.yticks(fontsize=8)
    plt.ylim(-bound,bound)
    ax.set_yticks(range(-int(bound),int(bound)))
    
    if normToGen==True:
        if normToMean==True:
            plt.ylabel("Mean & Gen-Norm'd Fitness Score (Temp vs. 28"+degree_sign+'C)', fontsize=fontsize-5)
        else:
            plt.ylabel('Generation-Normalized Fitness Score (Temp vs. 28'+degree_sign+'C)', fontsize=fontsize-5)
    else:
        if normToMean==True:
            plt.ylabel('Mean-Normalized Fitness Score (Temp vs. 28'+degree_sign+'C)', fontsize=fontsize)
        else:
            plt.ylabel('Raw Fitness Score (Temp vs. 28'+degree_sign+'C)', fontsize=fontsize)
    plt.xlabel('GOI', fontsize=fontsize)


    print('done plotting bar plot: '+figName)
    plt.tight_layout()

    plt.savefig(figName)
    
    return      



###RUN


sp_dict_37,sc_dict_37=ParseFitness(data_37)
sp_T_dict_37,sc_T_dict_37=calcTStats(sp_dict_37,sc_dict_37)
sp_M_dict_37,sc_M_dict_37=calcMean(sp_dict_37,sc_dict_37)

sp_norm_dict_37,sc_norm_dict_37=normToGen(sp_dict_37,gen_37),normToGen(sc_dict_37,gen_37)
sp_norm_T_dict_37,sc_norm_T_dict_37=calcTStats(sp_norm_dict_37,sc_norm_dict_37)


sp_dict_36,sc_dict_36=ParseFitness(data_36)
sp_T_dict_36,sc_T_dict_36=calcTStats(sp_dict_36,sc_dict_36)
sp_M_dict_36,sc_M_dict_36=calcMean(sp_dict_36,sc_dict_36)

sp_norm_dict_36,sc_norm_dict_36=normToGen(sp_dict_36,gen_36),normToGen(sc_dict_36,gen_36)
sp_norm_T_dict_36,sc_norm_T_dict_36=calcTStats(sp_norm_dict_36,sc_norm_dict_36)


###for just 36 and 37
####scatterplot of mean/stdev at 36 vs 37
###whole genome
##PlotScatter37_36(scatterFigName,sp_T_dict_36,sc_T_dict_36,sp_T_dict_37,sc_T_dict_37)
PlotScatter37_36("justMean_"+scatterFigName,sp_M_dict_36,sc_M_dict_36,sp_M_dict_37,sc_M_dict_37,justMean=True)
##PlotScatter37_36("normToGen_"+scatterFigName,sp_norm_T_dict_36,sc_norm_T_dict_36,sp_norm_T_dict_37,sc_norm_T_dict_37)
###only GOI
#PlotScatter37_36("goi_"+scatterFigName,sp_T_dict_36,sc_T_dict_36,sp_T_dict_37,sc_T_dict_37,onlyGOI=True)
PlotScatter37_36("goi_justMean_"+scatterFigName,sp_M_dict_36,sc_M_dict_36,sp_M_dict_37,sc_M_dict_37,justMean=True,onlyGOI=True)
#PlotScatter37_36("goi_normToGen_"+scatterFigName,sp_norm_T_dict_36,sc_norm_T_dict_36,sp_norm_T_dict_37,sc_norm_T_dict_37,onlyGOI=True)

##barplots: unnormalized
##Plot2Bars("sp_"+barFigName,sp_dict_36,sp_dict_37)
##Plot2Bars("sc_"+barFigName,sc_dict_36,sc_dict_37)
##
####barplots: fitness/generations
##Plot2Bars("sp_normToGen_"+barFigName,sp_norm_dict_36,sp_norm_dict_37,normToGen=True)
##Plot2Bars("sc_normToGen"+barFigName,sc_norm_dict_36,sc_norm_dict_37,normToGen=True)


#for all 3

sp_dict_38,sc_dict_38=ParseFitness(data_38)
sp_T_dict_38,sc_T_dict_38=calcTStats(sp_dict_38,sc_dict_38)
sp_M_dict_38,sc_M_dict_38=calcMean(sp_dict_38,sc_dict_38)

sp_norm_dict_38,sc_norm_dict_38=normToGen(sp_dict_38,gen_38),normToGen(sc_dict_38,gen_38)
sp_norm_T_dict_38,sc_norm_T_dict_38=calcTStats(sp_norm_dict_38,sc_norm_dict_38)

####barplots: unnormalized
##Plot3Bars("sp_"+barFigName3,sp_dict_36,sp_dict_37,sp_dict_38)
##Plot3Bars("sc_"+barFigName3,sc_dict_36,sc_dict_37,sc_dict_38)

####barplots: fitness/generations
#Plot3Bars("sp_normToGen_"+barFigName3,sp_norm_dict_36,sp_norm_dict_37,sp_norm_dict_38,normToGen=True)
#Plot3Bars("sc_normToGen_"+barFigName3,sc_norm_dict_36,sc_norm_dict_37,sc_norm_dict_38,normToGen=True)

####barplots: normalize fitness to mean (this did not work)
##Plot3Bars("sp_normToMean_"+barFigName3,normToMean(sp_norm_dict_36),normToMean(sp_norm_dict_37),normToMean(sp_norm_dict_38),normToGen=False,normToMean=True)
##Plot3Bars("sc_normToMean_"+barFigName3,normToMean(sc_norm_dict_36),normToMean(sc_norm_dict_37),normToMean(sc_norm_dict_38),normToGen=False,normToMean=True)
##
##Plot3Bars("sp_normToMean_normToGen_"+barFigName3,normToMean(sp_norm_dict_36),normToMean(sp_norm_dict_37),normToMean(sp_norm_dict_38),normToGen=True,normToMean=True)
##Plot3Bars("sc_normToMean_normToGen_"+barFigName3,normToMean(sc_norm_dict_36),normToMean(sc_norm_dict_37),normToMean(sc_norm_dict_38),normToGen=True,normToMean=True)

######scatterplot of mean/stdev at 37 vs 38
###whole genome
##PlotScatter37_38(scatterFigName_for38,sp_T_dict_37,sc_T_dict_37,sp_T_dict_38,sc_T_dict_38)
PlotScatter37_38("justMean_"+scatterFigName_for38,sp_M_dict_37,sc_M_dict_37,sp_M_dict_38,sc_M_dict_38,justMean=True)
##PlotScatter37_38("normToGen_"+scatterFigName_for38,sp_norm_T_dict_37,sc_norm_T_dict_37,sp_norm_T_dict_38,sc_norm_T_dict_38)
###only GOI
##PlotScatter37_38("goi_"+scatterFigName_for38,sp_T_dict_37,sc_T_dict_37,sp_T_dict_38,sc_T_dict_38,onlyGOI=True)
PlotScatter37_38("goi_justMean_"+scatterFigName_for38,sp_M_dict_37,sc_M_dict_37,sp_M_dict_38,sc_M_dict_38,justMean=True,onlyGOI=True)
##PlotScatter37_38("goi_normToGen_"+scatterFigName_for38,sp_norm_T_dict_37,sc_norm_T_dict_37,sp_norm_T_dict_38,sc_norm_T_dict_38,onlyGOI=True)
