import sys
import os
import numpy as np

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

'''
call as: python plot_genomecov.py d1.bedgraph d2.bedgraph[d1.bedgraph...]

input: bedgraph -d alignment(s) to reference genome
output: aligned plots
'''

bedgraphs = sys.argv[1:]

def ParseFromFiles(bedgraphs):
    '''
    Input: list of -d bedgraphs 
    Output: dict of {'begraph':[[positions],[coverages]]} dictionaries, sp1-16
    '''
    
    #build empty dictionary of chromosomes
    chrom_dict={}
##    bed_blanks={}
##    for d_bedgraph in bedgraphs:
##        bed_blanks[d_bedgraph]=[],[]
    for i in xrange(0,16):
        chrom_dict["sp"+str(i+1)]={}

    #build a coverage dictionary for each test file and add it to the chromosome dictionary  
    for d_bedgraph in bedgraphs:
        coverage_dict=ParseFromFile(d_bedgraph)
        for chromosome in chrom_dict:
            chrom_dict[chromosome][d_bedgraph]=coverage_dict[chromosome]
    print coverage_dict
    return chrom_dict

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

def plotAllCov(chromosome_dict, bedName):
    '''
    Input: a dictionary with bedgraphs as keys and pos,cov as values for a given chrom
    Output: saves a file per chromosome, with coverage plotted'''

    figName=bedName+".pdf"

    fig=plt.figure(figsize=(20,15))
    pl = plt.subplot(111)

    for d_bedgraph in chromosome_dict:
        chr_x=np.array(chromosome_dict[d_bedgraph][0])
        chr_y=np.array(chromosome_dict[d_bedgraph][1])
        pl.plot(chr_x,chr_y,label=d_bedgraph)
    
    pl.set_xlabel('position',fontsize=50)
    pl.set_ylabel('coverage',fontsize=50)

    plt.legend()
    plt.title(bedName)
    plt.savefig(figName, bbox_inches='tight', format='pdf', dpi=1000)

    return None


##RUN###

#make a directory for the image files

catNames=bedgraphs[0].split('.')[0]
for d_bedgraph in bedgraphs[1:]:
    catNames+="_"+d_bedgraph.split('.')[0]    
savePath=catNames+'_coveragePlots'

os.mkdir(savePath)


#Parse
chrom_dict=ParseFromFiles(bedgraphs)


#Plot
for chromosome in chrom_dict:
    bedName=savePath+'/'+chromosome+'_'+catNames
    plotAllCov(chrom_dict[chromosome],bedName)

#End
