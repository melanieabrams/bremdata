import regex
import numpy as np
import sys
import subprocess as sp
import copy
import pandas as pd


# PARAMETERS #    #Version 7-27-20 MBA
reverse_complement=True #set to True if barseq reads are the reverse complement of the NGS that goes into the TN-Seq mapping.

# HELP #

if len(sys.argv) == 1:
    print("USAGE: python map_PB_reads.py fastq_file ini_tnseq.fastq_pooled_reads out_directory")
    exit()
    
# INPUT # 

read_file = sys.argv[1]  
fastq_filename = read_file.split("/")[-1]
ini_tnseq_file = sys.argv[2]
out_dir = sys.argv[3]


# BEGIN FUNCTIONS #     

def getReverseComplement(bc):
    reverse_dict={'A':'T','T':'A','G':'C','C':'G'}
    revcomp_bc=""
    for base in bc:
        revcomp_bc=reverse_dict[base]+revcomp_bc
    return revcomp_bc

def parse_file(ini_tnseq_file, sep = '\t'):
    '''reads initial tnseq file as pd dataframe'''
    df=pd.read_csv(ini_tnseq_file)
    return df
    
def Filter_Reads(fastq_file):  ## filter out reads that do not have the universal primers

    bardict = {}

    wf = open(out_dir+fastq_filename+'_barcoded_reads','w')  # outfile for the trucated reads that will be mapped

    #organization: ....flanking_left-->barcode-->flanking_right
     
    flanking_left_pattern = regex.compile('(?e)(GATGTCCACGAGGTCTCT){e<=2}')  # Universal primer 1.  searches for this pattern allowing 2 mismatchs for sequencing errors
    flanking_right_pattern = regex.compile('(?e)(CGTACGCTGCAGGTCGAC){e<=2}')  # Universal primer 2.  searches for this pattern allowing 2 mismatchs for sequencing errors
    # counts for summary stats #

    tot_reads = 0.0  
    reads_with_bc = 0
    line_count = 0
    total_bases = 0
    
    f = open(fastq_file) # the file with the reads
    
    for line in f:

        # parses fastq file # 

        line_count+=1
        if line_count % 4 == 1:
            header = line
        elif line_count % 4 == 2:
            read = line.strip()
        elif line_count % 4 == 0:
            qual = line
            tot_reads+=1

            ###look for barcoding primers in the read sequence
            fl_match_data = flanking_left_pattern.search(read)  # searches for the primer pattern
            if fl_match_data == None:
                continue  #skips the read if it doesn't have universal primer 1

            fr_match_data = flanking_right_pattern.search(read)  # searches for the other primer pattern
            if fr_match_data == None:
                continue  #skips the read if it doesn't have a bc pattern

            barcode = read[fl_match_data.end():fr_match_data.start()]

            if reverse_complement==True: #added this 7/27/20 to account for NGS of TnSeq and BarSeq getting different strands
                barcode=getReverseComplement(barcode) 
            
            #calculate stats
            total_bases+=len(read)
            reads_with_bc+=1  ## counts reads that have the bc

            #add to dictionary
            if barcode in bardict:
                bardict[barcode]+=1
            else:
                bardict[barcode]=1
            
            wf.writelines(">"+header)
            wf.writelines(barcode+"\n")
            
    wf.close()    

    wf = open(out_dir+fastq_filename+"_mapping_stats", 'w')  # the file that will contain summary stats of the mapping and additional analysis

    # write mapping statistics" 

    print("total_reads: ", str(tot_reads))
    wf.writelines("total_reads: "+str(tot_reads)+"\n")
    print("reads with bc: ", str(reads_with_bc), " ("+str(100*reads_with_bc/tot_reads)[0:4]+"%)")
    wf.writelines("reads with bc: "+str(reads_with_bc)+" ("+str(100*reads_with_bc/tot_reads)[0:4]+"%)\n")


    wf.close()

    tot_parsed_reads = reads_with_bc

    return bardict

def assoc_with_barcode(ini_tnseq_file, bardict):
    '''creates a new csv file with counts from the barseq fastq
    associated with the mapped insert from the initial barseq tnseq'''
    df = pd.read_csv(ini_tnseq_file, sep = '\t') # make dataframe of initial mapping
    assoc_df= df.loc[df['ID'].isin(bardict.keys())] # pull out rows corresponding to barcodes found in this sequencing
    for index, row in assoc_df.iterrows():
        row_barcode = row['ID']
        count = bardict[row_barcode]
        assoc_df.at[index,'n'] = count #associate the metadata for that insert with the count from barseq
    assoc_df.to_csv(str(out_dir)+str('/')+str(fastq_filename)+'.barseq_pooled_reads', sep = '\t', index = False) #save as a new file in .fastq_pooled_reads format  
    return assoc_df
        
#### START PROGRAM ####

bardict = Filter_Reads(read_file)  ## filters out reads that don't have tn sequence, writes the genomic portion of the remaining reads to a new file

assoc_with_barcode(ini_tnseq_file, bardict)
