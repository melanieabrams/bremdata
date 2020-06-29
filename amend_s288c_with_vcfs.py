import sys
import vcf
import numpy as np

###USAGE: python amend_s288c_with_vcfs.py vcf h12files###
    #ex/ python ~/scripts/amend_s288c_with_vcfs.py ~/data/1WineEuropean_1011genomes/merged_1WineEuropean.vcf /usr2/people/mabrams/data/1WineEuropean_1011genomes/split_chrom/shs/all_samples/w1200_j25/*.h12_h2h1

    #ex/ python ~/scripts/amend_s288c_with_vcfs.py ~/data/1WineEuropean_1011genomes/split_chrom/merged_1WineEuropean_chromosome1.vcf /usr2/people/mabrams/data/1WineEuropean_1011genomes/split_chrom/shs/all_samples/w1200_j25/*.h12_h2h1



###PARAMETERS###

window_extension=3000 #number of bp to pull from upstream and downstream of the ORFs

reference_genome='/usr2/people/mabrams/Amended_Genomes/S288C/S288C_reference_genome_R64-1-1_20110203/S288C_reference_sequence_R64-1-1_20110203.fsa'

gff='/usr2/people/mabrams/Amended_Genomes/S288C/saccharomyces_cerevisiae_R64-1-1_20110208_withoutFasta.gff'
roman_numerals_in_gff=True

###FUNCTIONS###

def parseRefGenome(reference_genome):
    ''' Parses the reference genome
        Input: fasta formatted genome
        Output: dictionary {chrom:string},stripping newline characters'''
    
    roman_to_numerals={
        'I':'chr01','II':'chr02','III':'chr03','IV':'chr04','V':'chr05',
        'VI':'chr06','VII':'chr07','VIII':'chr08','IX':'chr09','X':'chr10',
        'XI':'chr11','XII':'chr12','XIII':'chr13','XIV':'chr14','XV':'chr15',
        'XVI':'chr16'}

    reference_dict={}
    f=open(reference_genome)
    seq="" #initialize seq
    for line in f:
        if line[0]=='>':
            if seq!="": #add the previous chromosome and header to te dictionary, reinitialize
                reference_dict[chrom]=seq
                seq=""
            if 'chromosome=' in line:
                chrom_rom=line.split('chromosome=')[1].strip()
                chrom_rom=chrom_rom[:-1]
                chrom=roman_to_numerals[chrom_rom]
            elif 'mitochondrion' in line:
                chrom='chrMito'
        else:
            line=line.strip()
            seq+=line
            
    #print(reference_dict['chr07'][(682566-1):687458]) #verified with ESP1, remember to subtract one for start to account for zero indexing, automatically comes off the end for exclusive list
    f.close()
    return reference_dict


def ParseFromGFF(gfffile):
    '''
    Parses SGD features flat file
    Input: SGD_features.tab file
    Output: dict of {chrom:{gene:[start,stop]}}
    '''
    roman_to_numerals={
        'chrI':'chr01','chrII':'chr02','chrIII':'chr03','chrIV':'chr04','chrV':'chr05',
        'chrVI':'chr06','chrVII':'chr07','chrVIII':'chr08','chrIX':'chr09','chrX':'chr10',
        'chrXI':'chr11','chrXII':'chr12','chrXIII':'chr13','chrXIV':'chr14','chrXV':'chr15',
        'chrXVI':'chr16','chrMito':'chrMito','2-micron':'chr2u'}
    ann_dict={}
    
    f = open(gfffile)
    lines=[]
    for line in f:
        if line[0]!='#': #skip header rows
            row_data=line.split("\t")
            #print(row_data)
            if roman_numerals_in_gff==True: #convert roman numerals if needed
                chrom=roman_to_numerals[row_data[0]]
            else:
                chrom=row_data[0]
            start=int(row_data[3])
            stop=int(row_data[4])
            info=row_data[8].split(";")
            yName=info[0].split('=')[1]
            #print(yName)
            if yName[0]=='Y' and len(yName)>5:
                if chrom in ann_dict:
                    ann_dict[chrom][yName]=[start,stop]
                else:
                    ann_dict[chrom]={yName:[start,stop]}
    f.close()
    
    return ann_dict

def helper_parseVCF(ann_dict):
    '''return location dicts to help vcf parser look only at loc of interest: {chrom: [[window_start_1,stop_1][w_s_2,stop_2]]}
       and list of chromosomes formatted like vcf'''
    #dict of windows of interest
    loc_dict={}
    for chrom in ann_dict:
        vcf_chr='chromosome'+str(int(chrom[3:])) #get chrom formatted like vcf
        for gene in ann_dict[chrom]:

            orf_start=ann_dict[chrom][gene][0] #orfs from annotation file
            orf_stop=ann_dict[chrom][gene][1]

            window_start=orf_start-window_extension #extend the window with the extension set in the PARAMETERs
            window_stop=orf_stop+window_extension

            start_index=window_start-1 # correcting for 0-1 indexing and exclusive slices
            stop_index=window_stop
            
            if vcf_chr in loc_dict:
                loc_dict[vcf_chr].append([start_index,stop_index,gene])
            else:
                loc_dict[vcf_chr]=[[start_index,stop_index,gene]]



    return loc_dict

def parseVCF(vcf_file,ann_dict):
    '''return nested dictionary {gene:{pos:[alleles]}}}'''

    #build loc_dict with helper function
    
    loc_dict=helper_parseVCF(ann_dict)

    #now go through vcf to sort variant calls by sample and gene
    
    vcf_dict={}
    for chrom in ann_dict:
        for gene in ann_dict[chrom]:
            vcf_dict[gene]={}

    vcf_reader = vcf.Reader(open(vcf_file, 'r'))
    for record in vcf_reader:
        chrom=record.CHROM
        if chrom in loc_dict: #only look at positions for chroms of interest
            print(chrom)
            pos=record.POS
            print(pos)
            #now see if position is in window of interest
            for window in range(len(loc_dict[chrom])):
                window_start=loc_dict[chrom][window][0]
                window_stop=loc_dict[chrom][window][1]
                if window_start<=pos<=window_stop:

                    #if it is, figure out what gene that was
                    gene=loc_dict[chrom][window][2]
                    
                    #add sample records to dictionary
                    for call in record.samples:
                        sample=call.sample
                        allele=call.gt_bases
                        if call.gt_type!=None and allele!=None:
                            alleles=allele.split('/')
                            if sample in vcf_dict[gene]:
                                vcf_dict[gene][sample][pos]=alleles
                            else:
                                vcf_dict[gene][sample]={pos:alleles}
                               
        else:
            print(chrom)

    for gene in vcf_dict:
        outfilename=gene+'_variants_log.txt'
        with open(outfilename,'w') as wf:
            for sample in vcf_dict[gene]:
                wf.writelines("#sample: "+sample+'\n')
                for pos in vcf_dict[gene][sample]:
                    wf.writelines(str(pos)+'\t'+ ','.join(vcf_dict[gene][sample][pos])+'\n')
                
    
    return vcf_dict



def amendGenes(gene, chrom, reference_dict, ann_dict, vcf_dict):
    '''output: {strain:seq} dict for a geoi'''

    strain_seqs={}

    #get window start and stop for gene
    orf_start=ann_dict[chrom][gene][0] #orfs from annotation file
    orf_stop=ann_dict[chrom][gene][1]

    window_start=orf_start-window_extension #extend the window with the extension set in the PARAMETERs
    window_stop=orf_stop+window_extension
    
    #get S288C sequence for that window
    start_index=window_start -1 
    stop_index=window_stop
    s288c_seq=reference_dict[chrom][start_index:stop_index]
    
    #get strain seqs
    for strain in vcf_dict[gene]:
        strain_seq_0="" #initialize the sequence
        strain_seq_1=""
        for i in range(len(s288c_seq)):
            pos=i+window_start
            if pos in vcf_dict[gene][strain]: 
                strain_seq_0+=vcf_dict[gene][strain][pos][0] #replace with variants
                strain_seq_1+=vcf_dict[gene][strain][pos][1]
            else:
                strain_seq_0+=s288c_seq[i]
                strain_seq_1+=s288c_seq[i]

        strain_seqs[strain+'_0']=strain_seq_0
        strain_seqs[strain+'_1']=strain_seq_1

    
    return strain_seqs, s288c_seq


def writeFasta(gene, amended_seqs, ref_seq):
    '''writes the amended gene into fasta format'''
    outfile=gene+'_amended.fa'
    n=80 #number of characters per line in output
    
    with open(outfile,'w') as wf:
        #write ref seq at the top
        wf.writelines('>ref(s288c)\n')
        seq=ref_seq
        split_n_chars=[seq[i:i+n]+'\n' for i in range(0, len(seq), n)]
        wf.writelines(split_n_chars)
        #then write following strains
        for header in amended_seqs:
            wf.writelines('>'+header+'\n')
            seq=amended_seqs[header]
            split_n_chars=[seq[i:i+n]+'\n' for i in range(0, len(seq), n)]
            wf.writelines(split_n_chars)
    return
    

### RUN ###

vcf_file=sys.argv[1]
h12files=sys.argv[2:]
save_prefix=h12files[0].split('_chromosome')[0]

#build annotation dict
print('parsing annotation file...')
ann_dict=ParseFromGFF(gff)

#build reference dict
print('parsing reference genome...')
reference_dict=parseRefGenome(reference_genome)


#get gene vars

print('loading variants from vcf...')
vcf_dict=parseVCF(vcf_file, ann_dict)

#amend gene
print('amending genes...')
for chrom in ann_dict:
    print('...from '+chrom)
    for gene in ann_dict[chrom]:
        amended_seqs, ref_seq=amendGenes(gene, chrom, reference_dict, ann_dict, vcf_dict)
        writeFasta(gene, amended_seqs, ref_seq)
    
