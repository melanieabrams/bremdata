import numpy as np
import sys


# PARAMETERS #    
ref_genome = "/usr2/people/mabrams/Amended_Genomes/Z1/Z1_global_ref.fa"  # looking for temp in Z1 since that's where I'd expect these genes to be unstable
gff = '/usr2/people/carlyweiss/SparScerinfo/YS2+CBS432+plasmid_clean' ## GFF to use for annotations.  

genes_of_interest=['YLR397C','YGR098C','YMR168C','YKR054C','YDR180W','YHR023W','YGR198W','YHR166C','YCR042C']



# BEGIN FUNCTIONS # 



def get_gff_dict(gff):
    '''gets chrom, start,end, strand for CDS annotations'''
    gff_dict = {}
    f = open(gff)

    for line in f:

        line = line.split("\t")
        chrom = line[1]
        gene = line[0]
        start = int(line[4])
        end = int(line[5])
        strand = line[3]
        ann_type = line[2]

        # put the annotation information in a dictionary #                                                                                                                                    
        if ann_type=='CDS':
            if gene[:2]=='sp': #because pulling from combined sp/sc annotation file
                gene=gene[2:]
                gff_dict[gene] = [chrom, start, end, strand]

    #print(gff_dict)
        
    return gff_dict

def parse_fasta(ref_genome):
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
    #print(chrom_dict['sp1'])

def translate(seq): 
    aa_table = { 
        'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M', 
        'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T', 
        'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K', 
        'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',                  
        'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L', 
        'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P', 
        'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q', 
        'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R', 
        'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V', 
        'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A', 
        'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E', 
        'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G', 
        'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S', 
        'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L', 
        'TAC':'Y', 'TAT':'Y', 'TAA':'STOP', 'TAG':'STOP', 
        'TGC':'C', 'TGT':'C', 'TGA':'STOP', 'TGG':'W', 
    } 
    protein =""
    
    if len(seq)%3 == 0: 
        for i in range(0, len(seq), 3): 
            codon = seq[i:i + 3]
            protein+= aa_table[codon]
    else:
        print('DNA sequence length remainder'+str(len(seq)%3))

    if protein[-4:]=='STOP':
        return protein[:-4]
    else:
        print('did not end in stop')
        return protein
            
def build_fa_goi(genes_of_interest, gff_dict, chrom_dict,mode='DNA'):
    yNames=""
    for gene in genes_of_interest:
        yNames+=(gene+"_")

    nt_table={'A':'T','T':'A','C':'G','G':'C'}

    wf = open(yNames+mode+'.fa','w')
    for gene in genes_of_interest:
        print(gene)
        
        
        #get annotation info for gene
        chrom=gff_dict[gene][0]
        start=gff_dict[gene][1]
        end=gff_dict[gene][2]
        strand=gff_dict[gene][3]

        #get DNA seq from genome
        seq=chrom_dict[chrom][start-1:end]
        if strand=='-':
##            if gene=='YPL174C':
##                print('here')
##                print(seq)
            reverse=seq[::-1]
            complement=''
            for nt in reverse:
                complement+=nt_table[nt]
            seq=complement

##        if gene=='YPL174C':
##            print('YPL174Cseq'+seq)
##            print(gff_dict['YPL174C'])
            
        #add DNA or translated sequence to outfile
        wf.writelines('>'+gene+'\n')
        if mode=='DNA':
            wf.writelines(seq+'\n')
        elif mode=='protein' or mode=='Protein':
            translated_seq=translate(seq)
            wf.writelines(translated_seq+'\n')
        else:
            print('need to specify DNA or protein')
        print()

    wf.close()
    return None

###START###


chrom_dict=parse_fasta(ref_genome)
gff_dict=get_gff_dict(gff)
build_fa_goi(genes_of_interest, gff_dict, chrom_dict,mode='protein')
        
