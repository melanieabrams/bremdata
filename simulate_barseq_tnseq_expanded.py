import regex
import numpy as np
import sys
import subprocess as sp
import random



##To-do: add known multimappers?  but I do trust blat.
##fix duplicate step so outfiles have the same number of lines.
##right now, it calls correct num w/ and without barcode, but thinks there's 136 duplicates instead of 100.
##to troubleshoot, decrease the number of everything to ~10 so I can manually check, maybe save outfiles that are more human-readable to do so
##also will need to find a way to get it to add a specified number of barcodes to existing genome positions (probably can replicate what I did for non-random barcodes)


# HELP #

##if len(sys.argv) == 1:
##    print("USAGE: python3 simulate_barseq_tnseq")
##    exit()
    
# PARAMETERS/INPUT # 

#genome file 
ref_genome ='/usr2/people/mabrams/Amended_Genomes/concatenated/D1373_Z1.fa' #will be pulled from this

#-FASTQ parameters: what should the reads themselves look like

num_without_barcode = 0 #number of reads to make without any barcode
num_with_barcode = 1500 #number of reads to make with a barcode
total_reads=num_without_barcode+num_with_barcode

bc_length = 20 # number of nt of bc
genome_nt = 50 #number of nt of genome in simulated read

#num_multimappers = 100 #number of barcodes that should DELIBERATELY be hard to distinguish between alleles - not implemented yet, maybe draw from a list of existing multimappers
num_nonunique_bc = 2 #number of barcodes that DELIBERATELY insert at multiple positions in the genome (can be higher by chance)
num_unique_bc=num_with_barcode-num_nonunique_bc
num_pos_multiple_barcodes=500 # number of genomic positions that DELIBERATELY have more than one distinct insertion.

#-Library parameters: how many distinct sets of insertions should there be, how much should they overlap
num_lib=2
num_dup_btn_tfn = 0 #number of reads that should DELIBERATELY have the same barcode at the same position in 2 the 2 libraries



# BEGIN FUNCTIONS #

def parse_fasta(ref_genome):
    '''returns a dictionary with each chromosome and the chromosome's sequence'''
    chrom_dict={}
    f=open(ref_genome)

    seq=""
    last_chrom=None
    for line in f:
        if line[0]== '>':
            chrom=line[1:].strip('\n')
            chrom_dict[chrom]=""
            if seq!="":
                chrom_dict[last_chrom]=seq
            seq=""
            last_chrom=chrom
                
        else:
            line=line.strip('\n')
            seq+=line

    return chrom_dict

def generate_bc(length=20):
    '''returns a random barcode of specified length'''
    random_bc = ''
    nucleotides = ['A','T','G','C']
    for i in range(length):
        random_bc+=random.choice(nucleotides)
    return random_bc

def get_genome_seq(genome_dict, length=50):
    '''pulls a random part of the genome'''
    chrom = random.choice(list(genome_dict.keys()))
    start = random.randint(0,len(genome_dict[chrom][:-abs(length)]))
    genome_seq = genome_dict[chrom][start:start+length]
    return genome_seq
    

def generate_read(genome_dict, barcode ='random', genome_seq = 'random'):
    '''returns a simulated barseq read with P5 and P7 adaptors (and Rd1 and Rd2 universal sequence primer) and a chunk of genome'''
    flanking_bc_left = 'AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCTNNNNNNAGTATGTAACCCTGATGTCCACGAGGTCTCT'
    #barcode can either be specified in input, random, or none
    #genome_seq can either be specified in input or random
    if barcode == 'none':
        if genome_seq == 'random':
            genome_seq = get_genome_seq(genome_dict, length=150)
        read = genome_seq
        return read, genome_seq, barcode
    else:
        if barcode == 'random':
            barcode= generate_bc()
        if genome_seq == 'random':
            genome_seq = get_genome_seq(genome_dict, length=genome_nt)
        flanking_bc_right = 'CGTACGCTGCAGGTCGACAACGTAAAACACATGCGTCAATTTTACGCATGATTATCTTTAACGTACGTCACAATATGATTATCTTTCTAGGGTTAA'
        after_genomic = 'AGATCGGAAGAGCACACGTCTGAACTCCAGTCACATCACGATCTCGTATGCCGTCTTC'
        read = flanking_bc_left + barcode + flanking_bc_right + genome_seq + after_genomic
        return read, genome_seq, barcode
        

def sim_lib(genome_dict, lib_num, dup=False, last_lib=['list with contents 0=read_seqs, 1=barcodes, 2=genome_seqs']):
    '''output simulated bar-seq data file'''
    wf = open('debug_simulatedn20_'+str(lib_num)+'.fq','w')  # outfile for the fake-barcoded reads that will be mapped

    read_seqs=[]
    genome_seqs=[]
    genomic_seqs_with_single_inserts=[]
    barcodes=[] #all barcodes
    distinct_barcodes=[] #non-redundant list of the barcodes

    if dup==False:
        num_this_tfn_barcode=num_with_barcode
        num_this_tfn_unique=num_unique_bc
    else: #decrease the number of new reads to generate by the number that you'll need to duplicate
        num_this_tfn_barcode=num_with_barcode - num_dup_btn_tfn
        num_this_tfn_unique=num_unique_bc - num_dup_btn_tfn

    for i in range(num_without_barcode):
        read, genome_seq, barcode = generate_read(genome_dict, barcode = 'none')
        read_seqs.append(read)
        genome_seqs.append(genome_seq)
        barcodes.append(barcode)

    for i in range(num_this_tfn_barcode):
        num_this_tfn_multiple_barcodes=0
        if i < num_this_tfn_unique:
            if num_this_tfn_multiple_barcodes<num_pos_multiple_barcodes and genomic_seqs_with_single_inserts!=[]:
                #add a randomly barcoded insert to an existing genomic position if conditions are met
                pos=random.choice(genomic_seqs_with_single_inserts)
                genomic_seqs_with_single_inserts.remove(pos)
                read, genome_seq, barcode = generate_read(genome_dict, barcode = 'random', genome_seq=pos)
            else:
                #otherwise, add a randomly barcoded insert to a new genomic position
                read, genome_seq, barcode = generate_read(genome_dict, barcode = 'random')
                genomic_seqs_with_single_inserts.append(genome_seq)
            read_seqs.append(read)
            genome_seqs.append(genome_seq)
            if barcode not in barcodes:
                distinct_barcodes.append(barcode)
            barcodes.append(barcode)
        else:
            #add specified number of repeat barcodes
            repeat_bc = random.choice(distinct_barcodes)
            read, genome_seq, barcode = generate_read(genome_dict, barcode = repeat_bc)
            read_seqs.append(read)
            genome_seqs.append(genome_seq)
            barcodes.append(barcode)
            
        
    if dup==True: #add entries from the last library if supposed to be duplicating
        for i in range(num_dup_btn_tfn):
            read_seqs.append(last_lib[0][i])
            barcodes.append(last_lib[1][i])
            genome_seqs.append(last_lib[2][i])
                    
    #write in an outfile formatted like the actual bb_merged barseq data (although shorter than those, since no extra stuff).
    read_count = 0
    for i in range(len(read_seqs)):
        read_count+=1
        wf.writelines("@simulateBarSeq_"+str(read_count)+'\n')
        wf.writelines(read_seqs[i]+'\n')
        wf.writelines('+'+'\n')
        wf.writelines('JJJ_LINE4_JJJ\n')
        
    wf.close()

    wf =open('debug_reads_simulated_'+str(lib_num)+'.fq','w')
    for read in read_seqs:
        wf.writelines(read+'\n')
    wf.close()
    
    wf =open('debug_barcodes_simulated_'+str(lib_num)+'.fq','w')
    for barcode in barcodes:
        wf.writelines(barcode+'\n')
    wf.close()
    
    wf =open('debug_genomeseqs_simulated_'+str(lib_num)+'.fq','w')
    for genomeseq in genome_seqs:
        wf.writelines(genomeseq+'\n')
    wf.close()

    return([read_seqs, barcodes, genome_seqs])

                  
#### START PROGRAM ####


genome_dict=parse_fasta(ref_genome)

for lib in range(num_lib):
    if lib==0:
        next_lib=sim_lib(genome_dict, lib)
    else:
        next_lib=sim_lib(genome_dict, lib, dup=True,last_lib=next_lib)
