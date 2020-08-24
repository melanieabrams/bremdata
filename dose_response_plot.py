import numpy as np
#import seaborn as sns
import matplotlib.pyplot as plt
csfont = {'fontname':'Arial'}

#PARAMETERS:

plot=True #set to True to plot the data
data_file='ESP1_swaps_raw_data.csv'
significance_cutoff=0.05
temps=[37.0,37.5,38.1,38.6,39.2,39.7,40.3,40.8]
dose_response_title='ESP1'
degree_sign='\u00b0'


OD_dict={
         'Sp,WT *':[0.367184856,0.259217943,0.146180476,0.063365507,0.02786568,0.017646095,0.022660294,0.014053092],
         'Sc,WT':[1.0,0.995205296,0.996371982,0.988869743,0.942526052,0.734575297,0.235043523,0.018382813],
         'ESP1 full swap *':[0.999896356,1.022675311,0.968660206,0.575370387,0.059984869,0.026462763,0.00883367,0.00089827]
        }


STD_dict={
         'Sp,WT *':[0.001028928,0.087060393,0.0338301,0.001307504,0.016115818,0.020159893,0.039547954,0.014222288],
         'Sc,WT':[0.0, 0.015381111, 0.039227154, 0.019300505, 0.002395222, 0.139807677, 0.216437943, 0.010276959],
         'ESP1 full swap *':[0.0, 0.015381111, 0.039227154, 0.019300505, 0.002395222, 0.139807677, 0.216437943, 0.010276959]
        }

###FNs###

def plotDoseResponse(temps, OD_dict, STD_dict):

    colors = [(1, 0.5, 0, 1),(0, 0.5, 1, 1), (0.05, 0.8, 0.1, 1),(0.05, 0.8, 0.1, 1), (0.05, 0.8, 0.1, 1),(0.05, 0.8, 0.1, 1), (0.05, 0.8, 0.1, 1)]
    fontsize=16
    strains=list(OD_dict.keys())
    x=np.array(temps) #x: num strains
    labels= [ ] #label: background, swap
    fig, ax=plt.subplots()

    ax.set(xlim=(37.0, 41.0), ylim=(0.0, 1.2))

    
    for i in range(len(strains)):
        strain=strains[i]
        labels.append(strain)
        ax.errorbar(x, np.array(OD_dict[strain]), yerr=STD_dict[strain],color=colors[i], marker='o',linewidth=1, markersize=6)

    xlabel='Temperature ('+degree_sign+'C)'
    
    #plt.title(dose_response_title,fontsize=fontsize,loc='center')
    plt.ylabel('Normalized Efficiency', fontsize=fontsize)
    plt.xlabel(xlabel, fontsize=fontsize)
    plt.legend(labels)

    plt.tight_layout()

    plt.savefig(dose_response_title+"_dose_response.png")    

    
    return                
                

        
        

        
###START###
print('done loading libraries')
plotDoseResponse(temps, OD_dict, STD_dict)
print('done plotting')
