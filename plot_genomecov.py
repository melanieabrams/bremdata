import sys
import os
import numpy as np

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

'''
call as: python plot_genomecov.py d1.bedgraph [d1.bedgraph...]

input: bedgraph -d alignment(s) to reference genome
output: aligned plots
'''

bedgraphs = sys.argv[1:]



#Test data: 

##sp1	1	0
##sp1	2	0
##sp1	3	0
##sp1	4	0
##sp1	5	0
##sp1	6	100
##sp1	7	0
##sp1	8	0
##sp1	9	0
##sp1	10	0



def ParseFromFile(d_bedgraph):
    '''
    Input: -d bedgraph (i.e., with a value per
    Output: {'chrom':[[positions],[coverages]]
    '''
    #build empty dictionary of chromosomes
    coverage_dict={}
    for i in xrange(0,16):
        coverage_dict["sp"+str(i+1)]=[],[]
         
    #add position and coverage for each base to a dictionary for that chromosome
    f = open(d_bedgraph)
    for line in f:
        row_data = line.strip().split("\t")
        chromosome = row_data[0]
        pos = int(row_data[1])
        cov = int(row_data[2])
        coverage_dict[chromosome][0].append(pos)
        coverage_dict[chromosome][1].append(cov)
    return coverage_dict

def plotCov(chromosome_vals, bedName):

    figName=bedName+".pdf"

    fig=plt.figure(figsize=(20,15))
    pl = plt.subplot(111)

    chr_x=np.array(chromosome_vals[0])
    chr_y=np.array(chromosome_vals[1])

    pl.plot(chr_x,chr_y)
    
    pl.set_xlabel('position',fontsize=50)
    pl.set_ylabel('coverage',fontsize=50)

    plt.title(bedName)
    plt.savefig(figName, bbox_inches='tight', format='pdf', dpi=1000)

    return None


##RUN###

for d_bedgraph in bedgraphs:
    #make a directory for the image files
    savePath=d_bedgraph+'coveragePlots'
    os.mkdir(savePath)
    
    #Parse
    coverage_dict=ParseFromFile(d_bedgraph)

    #Plot
    for chromosome in coverage_dict:
       bedName =savePath+'/'+d_bedgraph.split('.')[0]+"_"+chromosome
       plotCov(coverage_dict[chromosome],bedName)

