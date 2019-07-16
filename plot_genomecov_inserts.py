import sys
import os
import numpy as np
import matplotlib.pyplot as plt

'''
call as: python plot_genomecov.py file.fastq_pooled_reads

input: bedgraph -d alignment(s) to reference genome
output: aligned plots
'''

pooled_read_files = sys.argv[1:]
plot_all_on_one = True #set to false if you want it to plot to separate files

def ParseFromFile(pooled_reads):
    '''
    Input: _pooled_reads outfile from rh-seq pipeline
     tsv of: ID,scaffold,strand,location,annotation,n,rel_loc rel_prop,gene_length
    Output: {'chrom':[[positions],[coverages]]
    '''
    coverage_dict={}
        
    #add position and coverage for each base to a dictionary for that chromosome
    f = open(pooled_reads)
    next(f)
    for line in f:
        row_data = line.strip().split("\t")
        chromosome = row_data[1]
        if chromosome not in coverage_dict:
            coverage_dict[chromosome] = [],[]
        pos = int(row_data[3])
        cov = int(row_data[5])
        coverage_dict[chromosome][0].append(pos)
        coverage_dict[chromosome][1].append(cov)
    return coverage_dict

def plotCov(chromosome_vals, saveName):

    figName=saveName+".pdf"

    fig=plt.figure(figsize=(20,15))
    pl = plt.subplot(111)

    for i in range(len(chromosome_vals)):
        chr_x=np.array(chromosome_vals[i][0])
        chr_y=np.array(chromosome_vals[i][1])
        poolFile = pooled_read_files[i]
        pl.scatter(chr_x,chr_y,label=poolFile)
    
    pl.set_xlabel('position',fontsize=50)
    pl.set_ylabel('coverage',fontsize=50)

    plt.legend()
    plt.title(saveName)
    plt.savefig(figName, bbox_inches='tight', format='pdf', dpi=1000)

    plt.clf()
    return None


##RUN###

if plot_all_on_one == False:
    for pooled_read_file in pooled_read_files:
        #make a directory for the image files
        savePath=pooled_read_file+'coveragePlots'
        if not os.path.isdir(savePath):
            os.mkdir(savePath)
        
        #Parse
        coverage_dict=ParseFromFile(pooled_read_file)

        #Plot
        for chromosome in coverage_dict:
           saveName =savePath+'/'+pooled_read_file.split('.')[0]+"_"+chromosome
           plotCov([coverage_dict[chromosome]],saveName)

else:
    print(pooled_read_files)
    combined_name = pooled_read_files[0]
    for pooled_read_file in pooled_read_files[1:]:
        combined_name+='_'+pooled_read_file
    print(combined_name)
 
    #make a directory for the image files
    savePath=combined_name+'coveragePlots'
    if not os.path.isdir(savePath):
        os.mkdir(savePath)

    coverage_dicts = []
    for pooled_read_file in pooled_read_files:
        #Parse
        coverage_dict=ParseFromFile(pooled_read_file)
        coverage_dicts.append(coverage_dict)

    #Plot
    for chromosome in coverage_dicts[0]:
       saveName =savePath+'/'+pooled_read_file.split('.')[0]+"_"+chromosome
       chrom_dicts = []
       for i in coverage_dicts:
           chrom_dicts.append(i[chromosome])
       plotCov(chrom_dicts,saveName)

