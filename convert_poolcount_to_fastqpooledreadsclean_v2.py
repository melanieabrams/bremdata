import sys
from sys import argv
import pandas as pd
import numpy as np
import re
import string

### USAGE ###
# python organize_and_filter_genes.py output_folder filtered_insert_ratio_file.poolcounts 
#adds rel_loc 
#splits tech reps and converts poolcount file formatted for this into .fastq_pooled_read_clean for downstream rh-seq pipeline

##PARAMETERS##
gff='/usr2/people/mabrams/carlygithub/rh-seq/YS2+CBS433+plasmid_clean'
control_temp='28C' #assumes temp in name, e.g. 'Biorep1_28C_66'
#test_temp='4C'
#test_temp='39C'
test_temp='38'

##FUNCTIONS##
def parse_gff(gff):
    gff_dict = {}
    f = open(gff)

    for line in f:

        line = line.split("\t")
        chrom = line[1]
        gene = line[0]
        start = int(line[4])
        end = int(line[5])
        strand = line[3]
        type = line[2]

        # put the annotation information in a dictionary #                                                                                                                                    

        if chrom not in gff_dict:
            gff_dict[chrom] = {}

        gff_dict[chrom][gene] = [start, end, strand]
    f.close()
    return gff_dict


def parse_file(filename, sep='\t'):
    '''parses the poolcount file into a dataframe'''
    df = pd.read_csv(filename, sep=sep)
    file_string = str(filename)
    p = re.compile('[^/]+[^\.](?=\.)') # matches everything before the period and after the last slash, ie. the identifier of the file
    m = p.search(file_string) # finds regex in file name
    file_identifier = m.group() # prints match in string format
    return df, file_identifier

def retool_columns(df):
    #rename columns used in both but named differently
    df=df.rename(columns={'barcode':'ID','pos':'location','AlternateID':'annotation','CodingFraction':'rel_prop'})
    #get which columns have br n data
    br_columns=df.columns[8:]
    #drop unused columns
    df=df.drop(columns=['LocalGCpercent','NearestGene','Annotation'])
    #placeholders for columns without data
    df['strand']='seePoolFile'
    df['rel_loc']='notCalculated'
    df['gene_length']='notCalculated'
    #reorder columns to match Carly's format
    cols = df.columns.tolist()
    cols=cols[:4]+cols[5:-3]+[cols[-2]]+[cols[4]]+[cols[-1]]
    print(cols)
    #col_order=['ID','scaffold','location','annotation']+br_columns+['rel_loc','rel_prop','gene_length']
    df=df[cols] 
    #return reordered dataframe and the br columns
    return df,br_columns


def build_rename_dict(df,br_columns):
    '''builds a dictionary of new names for the bioreps named in the scheme:
    28A1 - biorep A of control temp, tech rep 1
    39A1 - birep A of test temp, tech rep 1'''

    control_brs=[]
    test_brs=[]
    t0_brs=[]
    other_brs=[]

    for col in br_columns:
        if control_temp in col:
            control_brs.append(col)
            
        elif test_temp in col:
            test_brs.append(col)

        elif 'T0_' in col:
            t0_brs.append(col)
                
        else:
            print('could not assign control or test temperature to column'+str(col))
            other_brs.append(col)

        
    print('control brs')
    print(control_brs)
    print('test brs')
    print(test_brs)
    #get alphabet for biorep identifiers
    alphabet=string.ascii_uppercase
    alphabet=alphabet.replace('E','') #remove E b/c bioreps named E introduce a bug downstream
    #generate renamin dictionary
    rename_dict={}
    for i in range(len(control_brs)):
        rename_dict[control_brs[i]]='28'+alphabet[i]+'1.fastq_pooled_reads_clean'
    for i in range(len(test_brs)):
        rename_dict[test_brs[i]]='39'+alphabet[i]+'1.fastq_pooled_reads_clean'
    for i in range(len(t0_brs)):
        rename_dict[t0_brs[i]]=t0_brs[i]+'.fastq_pooled_reads_clean'
    for i in range(len(other_brs)):
        rename_dict[other_brs[i]]='no outfile'

    print(rename_dict)
    
    return rename_dict

def split_br(df, rename_dict):
    '''splits out the data into bioreps according to the new filenames for downstream processing'''
    for br in rename_dict:
        if rename_dict[br]!='no outfile':
            other_brs=list(rename_dict.keys()) #get other brs
            #print('all_outfiles')
            #print(other_brs)
            print('making file for br:')
            print(br)
            #print('all outfiles remove br')
            other_brs.remove(br)
            #print(other_brs)

            br_df=df.copy() #drop other from a copy of df
            for other_br in other_brs: 
                br_df.drop(other_br, 1, inplace=True)

            br_df.rename(columns = {br:'n'}, inplace = True) 
       
            br_rename=rename_dict[br]
            br_df.to_csv(br_rename, sep='\t', index=False)
            
    return
        

##START
if __name__ == '__main__':
    print('...done importing libraries...')
    output_folder = argv[1].strip('/')
    poolcount_file = argv[2]
    print('...parsing gff...')
    gff_dict=parse_gff(gff)
    print('...parsing poolcounts file...')
    df,file_identifier=parse_file(poolcount_file)

    print('...renaming and adding columns...')
    df,br_columns=retool_columns(df)
    print(df.columns)
    print('...splitting by biorep...')
    rename_dict=build_rename_dict(df,br_columns)
    split_br(df,rename_dict)
    print('...done...')
    
