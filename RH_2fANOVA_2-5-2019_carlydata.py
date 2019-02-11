from sys import argv
from os import path

import math
import numpy as np
import matplotlib.pyplot as plt

import csv

from scipy import stats
import pandas as pd
from patsy import dmatrices

import statsmodels.api as sm
from statsmodels.formula.api import ols
import statsmodels.stats.multitest as smm

#to use python3.7:  #source activate /Applications/anaconda3

##modified from Carly Weiss's RH_6-6-2016.py RHA script by Melanie Abrams###
##Note: T0 refers to temperature 0, not time 0 (only final time temperature points used
## usage: add the files to the appropriate temperature categories in the table at the top, and then run as: python RH.py outdirectory

#ADD FILE NAMES & TEMPERATURES USED HERE:


temperatures_used = [28.0, 39.0] #temperatures MUST be in same order as files'

pools_T0 = ['fastq_pooledSRR7638932.fastq_pooled_reads','fastq_pooledSRR7638933.fastq_pooled_reads','fastq_pooledSRR7638934.fastq_pooled_reads',
            'fastq_pooledSRR7638935.fastq_pooled_reads','fastq_pooledSRR7638944.fastq_pooled_reads','fastq_pooledSRR7638945.fastq_pooled_reads',
            'fastq_pooledSRR7638946.fastq_pooled_reads','fastq_pooledSRR7638947.fastq_pooled_reads','fastq_pooledSRR7638950.fastq_pooled_reads',
            'fastq_pooledSRR7638951.fastq_pooled_reads','fastq_pooledSRR7638958.fastq_pooled_reads','fastq_pooledSRR7638959.fastq_pooled_reads'
            ]

pools_T1 = ['fastq_pooledSRR7638938.fastq_pooled_reads','fastq_pooledSRR7638939.fastq_pooled_reads',
            'fastq_pooledSRR7638940.fastq_pooled_reads','fastq_pooledSRR7638941.fastq_pooled_reads','fastq_pooledSRR7638942.fastq_pooled_reads',
            'fastq_pooledSRR7638943.fastq_pooled_reads','fastq_pooledSRR7638952.fastq_pooled_reads','fastq_pooledSRR7638953.fastq_pooled_reads',
            'fastq_pooledSRR7638954.fastq_pooled_reads','fastq_pooledSRR7638955.fastq_pooled_reads','fastq_pooledSRR7638956.fastq_pooled_reads',
            'fastq_pooledSRR7638957.fastq_pooled_reads'
            ]

sgdflatfile="/usr2/people/mabrams/Amended_Genomes/SGD_features.tab"



def ParseSGDFlat(sgdflatfile):
    '''
    Parses SGD features flat file
    Input: SGD_features.tab file
    Output: 
    '''

    descr_dict={}
    
    f = open(sgdflatfile)
    for line in f:
        row_data =line.split("\t")
        yName=row_data[3]
        descr=row_data[15]
        descr_dict[yName]=descr
    return descr_dict

def get_files(files_used, sep='\t', cols=[]):
        '''
        input: list of lists(temp categories) of .fastq_pooled replicates within that category
                .fasq_pooled are from Carly's map_and_pool blat py script, with the 9 following columns:
                        1.ID
                        2.scaffold
                        3.strand(+/-)
                        4.location
                        5.annotation
                        6.n (number of pooled reads per insert) *importapnt*
                        7.rel_loc
                        8.rel_prop
                        9.gene_length
        output: list of lists(temp categories) of list of file lines
        '''
        files = []
        for tempCat in files_used:
                catList=[]
                for replicate in tempCat:
                        lines = []
                        with open(replicate) as f:
                                f.readline()  # remove header
                                for line in f:
                                        line_list = line.rstrip().split(sep)
                                        lines.append([line_list[i] for i in cols])
                                catList.append(lines)
                files.append(catList)
                
##        for tempCat in files:
##                for replicate in tempCat:
##                        print(replicate)
##                
##         print(len(files),len(files[0]),len(files[0][0]))

                        
        return files #list of lists, where outermost layer is a list of temperature categories, next layer is a list of replicates within that category,
                     # next layer is a list of inserts, and the final layer is the insert's properties in the same order as the columns



def filter_inserts(files):
        '''remove inserts mapping to the plasmid and to non-coding regions, return filtered_files in same structure as get_files output'''
        filtered_files = []
        for tempCat in files:
                catList=[]
                for replicate in tempCat:
                        filtered_file = []
                        total_inserts = 0
                        total_map_plasmid = 0
                        number_inserts_not_genes = 0
                        too_few_reads = 0
                        for each_line in replicate:
                                total_inserts += 1
                                if each_line[1] == 'plasmid':  #ignore inserts mapping to plasmid
                                        total_map_plasmid += 1
                                        continue
                                if each_line[4] == 'NC':  #ignore inserts not mapping to genes
                                        number_inserts_not_genes += 1
                                        continue
                                else: 
                                        filtered_file.append(each_line)
                        catList.append(filtered_file)
                filtered_files.append(catList)

##        for tempCat in filtered_files:
##                for replicate in tempCat:
##                        print(replicate)
        
        return filtered_files  #list of files, each is a list of lines (as lists) that pass the filtering

def calc_total_mapped_reads(files):
        '''calculates and returns the total number of mapped reads per file (including reads mapping to plasmid, noncoding etc)'''
        total_mapped_reads_both = []
        for tempCat in files:
                catList = [ ]
                for x in tempCat:
                        total_mapped_reads = 0
                        for y in x:
                                total_mapped_reads+=float(y[5])
                        catList.append(total_mapped_reads)
                total_mapped_reads_both.append(catList)
        #print(total_mapped_reads_both)
        return total_mapped_reads_both  # a list of the number of total mapped reads per file in the order they were input in files_used

def make_data_dict(filtered_files,temperatures):  # the main chunk of code that generates a big dictionary of everything Melanie needed for a two-factor ANOVA
        '''output: multilayered dictionary where data_dict[gene_name][temperature][species] is a list of [normalized_read_count,rel_prop,ID,repCounter] for each insert'''
        data_dict = {}
        tempCounter=0 #counter variable to tell the script which of the temperatures from the temperatures list to assign the files in tempCat.  Could have made a dictionary, but was trying to match Carly's data structure AMAP.
        for tempCat in filtered_files: # for each T0 and T1
                temperature=temperatures[tempCounter]
                repCounter=0 #counter variable to index through replicates for normalization
                for replicate in tempCat: #for each replicate 
                        for insert in replicate:  # for each insert
                                gene_name = insert[4].strip('sc').strip('sp')
                                species = insert[4][0:2]  # defines the species as either sp or sc, depending on the prefix in the gene name
                                ID = str(insert[0])
                                rel_prop = float(insert[7])

                                if gene_name not in data_dict:
                                        data_dict[gene_name] = {temperatures[0]:{'sc':[],'sp':[]},temperatures[1]:{'sc':[],'sp':[]}}
                                        
                                normalized_read_count = float(insert[5]) / total_mapped_reads_both[tempCounter][repCounter]
                                data_dict[gene_name][temperature][species].append([normalized_read_count, rel_prop, ID,repCounter])
                                
                        repCounter+=1
                tempCounter+=1


        #print(data_dict)
        return data_dict


def how_many_inserts_per_gene(data_dictionary):
        '''
        input: input: dictionary with layers  gene, temperature, species [number of reads, insertID, replicate]
        output: writes a file of total number and allele-specific # of inserts of the input data_dictionary
        '''

        wf = open(out_dir+"number_inserts_per_gene.txt", "w")
        wf.writelines("Gene name\tSp inserts\tSc inserts\tTotal inserts\n")

        for gene_name in data_dictionary:
                total_sc_inserts=0
                total_sp_inserts=0
                for temperature in data_dictionary[gene_name]:
                        for species in data_dictionary[gene_name][temperature]:
                                for insert in data_dictionary[gene_name][temperature][species]:
                                        if species == 'sc':
                                                total_sc_inserts+=1
                                        if species == 'sp':
                                                total_sp_inserts+=1
                total_number_inserts=total_sc_inserts+total_sp_inserts
                wf.writelines(str(gene_name)+"\t"+str(total_sp_inserts)+"\t"+str(total_sc_inserts)+"\t"+str(total_number_inserts)+"\n")

        wf.close()

def filter_data_dict(data_dictionary, number_inserts_per_allele_needed=5):
        '''
       input: input: dictionary with layers  gene, temperature, species [number of reads, insertID, replicate]
       output: [0]  dictionary in the same format, having filtered out genes that don't have enough unique insertions in BOTH the species' alleles of that gene across ALL temperatures/replicates
               note: does NOT filter for a minimum number of reads per insert -- can add that later if I need to, but those thresholds were 1 in the paper
               [1]  dictionary of log ratios of inserts
        '''

        filtered_data_dict = {}
        removed_inserts_position = {} # maybe one day I will script it so that the ones that I don't use are added to another dict. not now though ><

        counterDict = {} # {gene: [total_sc_inserts, total_sp_inserts, [rel_props (i.e., unique inserts)]]}
        ratio_dict= {} #{gene: total_sc_inserts_T0, total_sc_inserts_T1, total_sp_inserts_T0, total_sp_inserts_T1, sc_ratio, sp_ratio}

			
        for gene_name in data_dictionary:
                counterDict[gene_name]=[0,0,[]]
                ratio_dict[gene_name]=[0,0,0,0]
                tempCounter=0
                for temperature in data_dictionary[gene_name]:
                        for species in data_dictionary[gene_name][temperature]:
                                if species == 'sc':
                                        for insert in data_dictionary[gene_name][temperature][species]: # iterate over each insert
                                                if insert[1] not in counterDict[gene_name][2]: #check whether that exact insert has already shown up
                                                        counterDict[gene_name][0]+=1
                                                        counterDict[gene_name][2].append(insert[1])
                                                        if tempCounter==0:
                                                                ratio_dict[gene_name][1]+=1
                                                        if tempCounter==1:
                                                                ratio_dict[gene_name][0]+=1
                                                        
                                if species == 'sp':
                                        for insert in data_dictionary[gene_name][temperature][species]: # iterate over each insert
                                                if insert[1] not in counterDict[gene_name][2]: #check whether that exact insert has already shown up
                                                        counterDict[gene_name][1]+=1
                                                        counterDict[gene_name][2].append(insert[1])
                                                        if tempCounter==0:
                                                                ratio_dict[gene_name][3]+=1
                                                        if tempCounter==1:
                                                                ratio_dict[gene_name][2]+=1
                        tempCounter+=1
                #print(ratio_dict)
                try:
                        filler_number_T0 = 1/(ratio_dict[gene_name][0]+ratio_dict[gene_name][2]) #if no reads were found for that insert in either T0 or T1, need a "filler" number to use to rep 1 read so that not dividing over 0
                except ZeroDivisionError:
                        filler_number_T0=0.000000000001
                try:
                        filler_number_T1 = 1 /(ratio_dict[gene_name][1]+ratio_dict[gene_name][3])
                except ZeroDivisionError:
                        filler_number_T1=0.000000000001
                if ratio_dict[gene_name][0]==0:
                        ratio_dict[gene_name][0]=filler_number_T0
                if ratio_dict[gene_name][1]==0:
                        ratio_dict[gene_name][1]=filler_number_T1
                if ratio_dict[gene_name][2]==0:
                        ratio_dict[gene_name][2]=filler_number_T0
                if ratio_dict[gene_name][3]==0:
                        ratio_dict[gene_name][3]=filler_number_T1

                ratio_dict[gene_name].append(math.log(ratio_dict[gene_name][1]/ratio_dict[gene_name][0],2)) #adds the log ratios
                ratio_dict[gene_name].append(math.log(ratio_dict[gene_name][3]/ratio_dict[gene_name][2],2))


        #print(counterDict)        
        for gene_name in counterDict:
                if counterDict[gene_name][0]>=number_inserts_per_allele_needed and counterDict[gene_name][1]>= number_inserts_per_allele_needed:
                        filtered_data_dict[gene_name]=data_dictionary[gene_name]

        #print(counterDict['YMR130W'])
        #print(ratio_dict['YMR130W'])
        #print(filtered_data_dict)
        return filtered_data_dict, ratio_dict

def make_data_frame(filtered_data_dict,files_used, temperatures):
        '''
        input: dictionary with layers  gene, temperature, species [number of reads, insertID, replicate]
        output: dictionary {gene:pandas dataframe} in the following structure

        for a given gene

        (row_label)     FILE    TEMPERATURE     SPECIES         NUM READS
        28_1_sc         28_1    28              sc              n_28_1_sc
        28_1_sp         28_1    28              sp              n_28_1_sp
        28_2_sc         28_2    28              sc              n_28_2_sc
        28_2_sp         28_2    28              sp              n_28_2_sp
        28_3_sc         28_3    28              sc              n_28_3_sc
        28_3_sp         28_3    28              sp              n_28_3_sp
        39_1_sc         39_1    39              sc              n_39_1_sc
        39_1_sp         39_1    39              sp              n_39_1_sp
        39_2_sc         39_2    39              sc              n_39_2_sc
        39_2_sp         39_2    39              sp              n_39_2_sp
        39_3_sc         39_3    39              sc              n_39_3_sc
        39_3_sc         39_3    39              sp              n_39_3_sc

        use these dataframes to feed into the ols model-maker for two-factor ANOVA test
        '''
        all_dataframes={}
        
        #build empty dictionary containing ALL files and all temps
        dfdict={}
        for tempIndex in range(len(files_used)):
                temp = temperatures[tempIndex]
                repCounter=0
                for repFile in files_used[tempIndex]:
                        repCounter+=1
                        fileVal=str(int(temp))+'_'+str(repCounter)
                        rowVal_sc=str(int(temp))+'_'+str(repCounter)+'_'+'sc'
                        rowVal_sp=str(int(temp))+'_'+str(repCounter)+'_'+'sp'
                        dfdict[rowVal_sc]=[fileVal,temp,'sc',0.0]
                        dfdict[rowVal_sp]=[fileVal,temp,'sp',0.0]
        #print(dfdata)

        #make dataframe for the gene from a dictionary in the dfdict template, add it to all_dataframes                
        for gene_name in filtered_data_dict:
                dfgenedata=dfdict.copy()
                for temperature in filtered_data_dict[gene_name]:
                        for species in filtered_data_dict[gene_name][temperature]:
                                for insert in filtered_data_dict[gene_name][temperature][species]:
                                        replicate=1+insert[-1]
                                        num_reads=insert[0]
                                        this_row=str(int(temperature))+'_'+str(replicate)+'_'+species
                                        dfgenedata[this_row][3]+=num_reads
                                
                dfgene=pd.DataFrame.from_dict(dfgenedata,orient='index', columns=['FILE','TEMPERATURE','SPECIES','NREADS'])
                all_dataframes[gene_name]=dfgene
##
##        interesting_genes=['YOR326W','YOR375C','YOR363C', 'YML065W','YMR148W','YMR001C']
##        interesting_genes=['YGR098C']
##        for gene in interesting_genes:
##            saveName=out_dir+gene+'.txt'
##            all_dataframes[gene].to_csv(saveName,header=True,index=False,sep='\t',mode='a')
##
##                
##        for df in all_dataframes:
##                print(df)
##                print(filtered_data_dict[df])
##                print(all_dataframes[df])
##                print()
                
        return all_dataframes

class HemizygoteComparison(object):
        '''now doing the RH test'''
        FIELDNAMES = ["Gene", "Fvalue", "pvalue","padjusted", "Number Sc Inserts", "Number Sp Inserts","Aggregate Sc Fitness", "Aggregate Sp Fitness","SGD"]

        def __init__(self, gene_name, Fval, pval, n_sc_inserts, n_sp_inserts, agg_sc_fitness, agg_sp_fitness,sgd):
                self.gene_name = gene_name
                self.Fval = Fval
                self.pvalue = pval
                self.adjusted_pval = None
                self.number_sc_inserts = n_sc_inserts
                self.number_sp_inserts = n_sp_inserts
                self.aggregate_sc_fitness = agg_sc_fitness
                self.aggregate_sp_fitness = agg_sp_fitness
                self.sgd=sgd
 

        def as_dict(self):
                return {'Gene': self.gene_name, 
                                'Fvalue': self.Fval,
                                'pvalue': self.pvalue,
                                'padjusted': self.adjusted_pval
                                ,'Number Sc Inserts': self.number_sc_inserts,
                                'Number Sp Inserts': self.number_sp_inserts,
                                'Aggregate Sc Fitness': self.aggregate_sc_fitness,
                                'Aggregate Sp Fitness': self.aggregate_sp_fitness,
                                'SGD': self.sgd
                                }


def reciprocal_hemizygote_test(dataframes_dict,ratio_dict):
        '''executes RH test on each gene'''
        
        sgd_dict=ParseSGDFlat(sgdflatfile)
        
        RHlist = []
        for gene_name in dataframes_dict:
                df = dataframes_dict[gene_name]
                species_temp_lm = ols('NREADS ~ C(SPECIES)*C(TEMPERATURE)',df).fit()
                RH_test = sm.stats.anova_lm(species_temp_lm,typ=2)
                pvalue=RH_test['PR(>F)'][2] 
                Fval=RH_test['F'][2]
##                if gene_name=='YGR098C':
##                    print(gene_name)
##                    print(species_temp_lm)
##                    print(RH_test, pvalue,Fval)
        
        ##for gene_name in ratio_dict:
                number_sc_inserts = 0
                number_sp_inserts = 0
                if type(ratio_dict[gene_name][0])==int:
                        number_sc_inserts+=ratio_dict[gene_name][0]
                if type(ratio_dict[gene_name][1])==int:
                        number_sc_inserts+=ratio_dict[gene_name][1]
                if type(ratio_dict[gene_name][2])==int:
                        number_sp_inserts+=ratio_dict[gene_name][2]
                if type(ratio_dict[gene_name][3])==int:
                        number_sp_inserts+=ratio_dict[gene_name][3]
                aggregate_sc_fitness = ratio_dict[gene_name][4]
                aggregate_sp_fitness = ratio_dict[gene_name][5]
                comparison_data = HemizygoteComparison(gene_name, Fval, pvalue,
                        number_sp_inserts, number_sc_inserts,
                        aggregate_sp_fitness, aggregate_sc_fitness,sgd_dict[gene_name])

                RHlist.append(comparison_data)
	
        pvalue_list = [cd.pvalue for cd in RHlist]

        #stats2 = importr('stats') # for using rpy2
        #adjusted_pvalues = stats2.p_adjust(FloatVector(pvalue_list), method = 'BH')  #do multiple testing correction with BH method
        fdrbh_output = smm.multipletests(pvalue_list, method='fdr_bh')
        adjusted_pvalues = fdrbh_output[1].tolist()
        for cd, padjusted in zip(RHlist, adjusted_pvalues):
                cd.adjusted_pval = padjusted


        with open(out_dir+"testRHresults.txt", "w") as wf:
                writer = csv.DictWriter(wf, HemizygoteComparison.FIELDNAMES, dialect='excel-tab')
                writer.writeheader()
                for cd in RHlist:
                        writer.writerow(cd.as_dict())





####  Running the script  ####
print("libraries loaded, getting started..")
files_used = [pools_T0, pools_T1]
out_dir=argv[1]+'_'

#parse files
files = get_files(files_used, cols=range(0, 9))
print("files loaded...")
filtered_files=filter_inserts(files)
print("filtered out intergenic regions..")
total_mapped_reads_both = calc_total_mapped_reads(files)
print("wrote total_mapped_reads file...")
data_dict = make_data_dict(filtered_files,temperatures_used)
print("parsing and filtering data...")
filtered_data_dict, filtered_ratio_dict = filter_data_dict(data_dict, number_inserts_per_allele_needed=5) 
how_many_inserts_per_gene(filtered_data_dict) ## calculate how many inserts per gene in each allele that pass all the thresholds
dataframe_dict = make_data_frame(filtered_data_dict,files_used,temperatures_used)
print("ready to perform RH test")
reciprocal_hemizygote_test(dataframe_dict,filtered_ratio_dict)
print("done!")






