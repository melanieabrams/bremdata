import regex
import numpy as np
import sys
import subprocess as sp
import random


# HELP #

if len(sys.argv) == 1:
    print("USAGE: python3 simulate_barseq_tnseq out_directory fastq_file1 fastq_file2...")
    exit()
    
# INPUT # 

num_orig = 20 #number of unmodified reads to preserve.  This will make sure my modified version for barseq of map-and-blat can still filter out reads w/o Tn.
num_new = 1000 # number of new reads
num_duplicate = 100 #number of new reads with duplicate barcodes
bc_length = 20 # number of nt of bc
genome_nt = 50 #number of nt of genome in simulated read

# BEGIN FUNCTIONS #

def generate_bc(length=20):
    '''returns a random barcode of specified length'''
    random_bc = ''
    nucleotides = ['A','T','G','C']
    for i in range(length):
        random_bc+=random.choice(nucleotides)
    return random_bc

def generate_read(genome_seq,barcode ='random'):
    '''returns a simulated barseq read with P5 and P7 adaptors (and Rd1 and Rd2 universal sequence primer) and a chunk of genome'''
    flanking_bc_left = 'AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCTNNNNNNAGTATGTAACCCTGATGTCCACGAGGTCTCT'
    if barcode == 'random':
        barcode= generate_bc()
    flanking_bc_right = 'CGTACGCTGCAGGTCGACAACGTAAAACACATGCGTCAATTTTACGCATGATTATCTTTAACGTACGTCACAATATGATTATCTTTCTAGGGTTAA'
    after_genomic = 'AGATCGGAAGAGCACACGTCTGAACTCCAGTCACATCACGATCTCGTATGCCGTCTTC'
    read = flanking_bc_left + barcode + flanking_bc_right + genome_seq + after_genomic
    return read
    

def AddN20(fastq_file):  ## add a random N20 plus primers to reads with transposon, so that normal Tn-Seq data looks like it was made with a barcode

    wf = open(out_dir+fastq_filename+'_simulatedn20','w')  # outfile for the fake-barcoded reads that will be mapped
        
    line_count = 0
    tn_count = 0

    with open(fastq_file) as f:
        head = [next(f) for x in range(4*(num_orig+num_new))]
    
    for line in head:
        line_count +=1
        if line_count % 4 == 1:
            header = line 
        elif line_count % 4 == 2:
            read = line.strip()
        elif line_count % 4 == 0:
            nt_from_read =read[75:75+genome_nt]
            if line_count >4*num_orig:
                if line_count>4*(num_new-num_duplicate): 
                    read = generate_read(nt_from_read,barcode='random')
                else:
                    read = generate_read(nt_from_read,barcode='TATTGGAAAACTATAGGGAC')
            wf.writelines(">simulatedBarSeq"+header)
            wf.writelines(read+"\n")            

                  

#### START PROGRAM ####

out_dir = sys.argv[1]
read_files = sys.argv[2:]  

for read_file in read_files:
    fastq_filename = read_file.split("/")[-1]
    AddN20(read_file)
