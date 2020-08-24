import numpy as np
#import seaborn as sns
import matplotlib.pyplot as plt
csfont = {'fontname':'Arial'}

#PARAMETERS:

degree_sign='\u00b0'

#dict: {strain_temp: [large,med,small]}
lms_dict={'Sp,WT_28':[0,0.25,0.75],
          'Sp,WT_39':[0.611111111,0.368055556,0.020833333],
          'Sc,WT_28':[0.338509317,0.431677019,0.229813665],
          'Sc,WT_39':[0.598214286,0.125,0.276785714],
          'ESP1 full swap_28':[0.339215686,0.598039216,0.062745098],
          'ESP1 full swap_39':[1.0,0.0,0.0]}


def plotDyads(lms_dict):

    N=6
    ind=np.arange(N)
    
    width=0.6
    #colors = [(1, 0.5, 0, 1),(0, 0.5, 1, 1), (0.05, 0.8, 0.1, 1),(0.05, 0.8, 0.1, 1), (0.05, 0.8, 0.1, 1),(0.05, 0.8, 0.1, 1), (0.05, 0.8, 0.1, 1)]
    fontsize=16
    strains_temps=list(lms_dict.keys())
    
    #x=np.array(temps) 
    xlabels= [ ]

    fig, ax=plt.subplots(figsize=(10,6))

    top_side = ax.spines["top"]
    top_side.set_visible(False)

    #ax.set(xlim=(37.0, 41.0), ylim=(0.0, 1.2))

    lg=[]
    med=[]
    sm=[]
    for i in range(len(strains_temps)):
        #print(strains_temps[i].split('_'))
        strain,temp=strains_temps[i].split('_')
        xlabels.append(strain+'\n'+temp+degree_sign+'C')
        lg.append(lms_dict[strains_temps[i]][0])
        med.append(lms_dict[strains_temps[i]][1])
        sm.append(lms_dict[strains_temps[i]][2])

    
    lg_and_med=[]
    for i in range(len(lg)):
        lg_and_med.append(lg[i]+med[i])


    
    ax.bar(xlabels, lg, width, label='large',color='purple')
    ax.bar(xlabels,med,width,bottom=lg, label='medium',color='brown')
    ax.bar(xlabels,sm,width,bottom=lg_and_med, label='small')


    #ax.legend(labels)

    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    
    #plt.title('Dyad Barplot',fontsize=fontsize,loc='center')
    plt.ylabel('Fraction of dyads', fontsize=fontsize)
    plt.xlabel('Strain', fontsize=fontsize)
    

    plt.tight_layout()

    plt.savefig("dyad_barplot.png")    

    
    return                
                

        
        

        
###START###
print('done loading libraries')
plotDyads(lms_dict)
print('done plotting')
