import sys
import numpy as np
import matplotlib.pyplot as plt

'''
plot sweepfinder 2 outputs
names chromosome assuming file in form: Bergstrom2014_Spar_1.sf2  
'''

filenames = sys.argv[1:]
plot_all_on_one = True #set to false if you want it to plot to separate files

def ParseSF2(sf2_file):
    '''
    Input: outfile from sweepfinder 2
    Output: 
    '''
    pos_list=[]
    LR_list=[]
        
    #add position and coverage for each base to a dictionary for that chromosome
    f = open(sf2_file)
    next(f)
    for line in f:
        row_data = line.strip().split("\t")
        pos_list.append(float(row_data[0]))
        LR_list.append(float(row_data[1]))

    f.close()

    return pos_list, LR_list

def plotLR(pos_list, LR_list, prefix):

    fig=plt.figure(figsize=(20,15))
    pl = plt.subplot(111)

    x=np.array(pos_list)
    y=np.array(LR_list)
    pl.scatter(x,y)
    
    pl.set_xlabel('position',fontsize=50)
    pl.set_ylabel('SF2 LR',fontsize=50)

    plt.title(prefix+' LR plot')
    plt.savefig(figname, bbox_inches='tight', format='pdf', dpi=1000)

    plt.clf()
    return None


##RUN###

for sf2file in filenames:
    pos_list, LR_list=ParseSF2(sf2file)
    prefix=sf2file.split('.')[0]
    figname=prefix+'_LRplot.pdf'
    plotLR(pos_list, LR_list, prefix)
