import sys
import vcf
import numpy as np

#USAGE: python amend_s288c_with_vcfs.py vcf h12files
#ex/ python ~/scripts/amend_s288c_with_vcfs.py ~/data/1WineEuropean_1011genomes/merged_1WineEuropean.vcf *.h12_h2h1

#PARAMETERS
goi=['YLR397C','YGR098C', 'YMR168C','YKR054C','YDR180W',
     'YHR023W','YGR198W','YHR166C','YCR042C','YPL174C']

#goi=['YGR098C']

reference_genome='/usr2/people/mabrams/Amended_Genomes/S288C/S288C_reference_genome_R64-1-1_20110203/S288C_reference_sequence_R64-1-1_20110203.fsa'

gff='/usr2/people/mabrams/Amended_Genomes/S288C/saccharomyces_cerevisiae_R64-1-1_20110208_withoutFasta.gff'
roman_numerals_in_gff=True

def parseRefGenome(reference_genome):
    '''return dictionary {chrom:string},stripping newline characters'''
    
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

def ParseH12(h12_file):
    '''
    Input: outfile from selscan
    Output: a dictionary of winodws tested, in the form {center: right, left}
    '''
    pos_dict={}

    #add position and values for each base to a dictionary for that chromosome
    f = open(h12_file)
    next(f)
    for line in f:
        row_data = line.strip().split("\t")
        center_coord=int(row_data[0])
        left_coord=int(row_data[1])
        right_coord=int(row_data[2])
        G1=float(row_data[6])
        pos_dict[center_coord]=[left_coord, right_coord,G1]
        
        #columns:
            #1ctrcoord (index 0)
            #2leftcoord
            #3rightcoord
            #4K
            #55hapfreqspec
            #6strainnum
            #7H1 (index 6)
            #8H2 (index 7)
            #9H12 (index 8)
            #10H2/H1 (index 9)
            #11H123 (index 10)
 
    f.close()

    return pos_dict 


def getStartStop(ann_dict, chrom_dict):
    '''get start and stop for the G12 window around each gene, as measured by ranking script'''
    outfileName=save_prefix+'_gene_windows.txt'
    G1_dict={}
    gene_G1s=[]
    for chrom in ann_dict:
        pos_dict=chrom_dict[chrom] # centers of windows for G1
        for gene in ann_dict[chrom]:
            G1s=[]
            coords=[]
            start=ann_dict[chrom][gene][0] #annotated gene start
            stop=ann_dict[chrom][gene][1] #annotated gene stop
            for pos in pos_dict:
                if start<pos<stop:
##                    if gene=="YAL056W":
##                        print(gene)
##                        print("H12 list")
##                        print(H12_list[i])
##                        print("H2/H1 ratios")
##                        print(ratio_H2H1_list[i])
                    coords+=pos_dict[pos][:-1] #add the right and left coordinate to t
                    G1s.append(pos_dict[pos][-1])
                    

            if len(G1s)>0:
                gene_G1=np.mean(G1s)
                gene_G1s.append(gene_G1)
                window_start=min(coords)
                window_stop=max(coords)
                G1_dict[gene]=[gene_G1, chrom, start, stop, window_start, window_stop]


    with open(outfileName, 'w') as wf:
        wf.writelines('gene\taverageG1\tchrom\tstart\tstop\twindow_start\twindow_stop\n')
        for k in sorted(G1_dict, key=G1_dict.get, reverse=True):
            #print(gene)
            gene=k
            avgG1=str(G1_dict[gene][0])
            chrom=str(G1_dict[gene][1])
            start=str(G1_dict[gene][2])
            stop=str(G1_dict[gene][3])
            window_start=str(G1_dict[gene][4])
            window_stop=str(G1_dict[gene][5])
            wf.writelines(gene+'\t'+avgG1+'\t'+chrom+'\t'+start+'\t'+stop+'\t'+window_start+'\t'+window_stop+'\n')
            
    return G1_dict

def helper_parseVCF(coord_dict):
    '''return location dicts to help vcf parser look only at loc of interest: {chrom: [[window_start_1,stop_1][w_s_2,stop_2]]}
       and list of chromosomes formatted like goi'''
    #dict of windows of interest
    loc_dict={}
    for gene in goi:
        chrom=coord_dict[gene][1]
        chrnum=int(chrom[3:])
        chrom='chromosome'+str(chrnum) #formatted like vcf
        
        window_start=coord_dict[gene][4]
        window_stop=coord_dict[gene][5]

        start_index=window_start-1 # correcting for 0-1 indexing and exclusive slices
        stop_index=window_stop
        
        if chrom in loc_dict:
            loc_dict[chrom].append([start_index,stop_index,gene])
        else:
            loc_dict[chrom]=[[start_index,stop_index,gene]]



    return loc_dict

def parseVCF(vcf_file,coord_dict):
    '''return nested dictionary {gene:{pos:[alleles]}}}'''

##    if True: #for debugging
##        return {'YGR098C': {'AAA': {687458: ['G', 'T']}}}


    #helper

    
    loc_dict=helper_parseVCF(coord_dict)
    #print(loc_dict)

    #now    
    
    vcf_dict={}
    for gene in goi:
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
                                #print(vcf_dict)
                                #return vcf_dict #for troubleshooting
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



def amendGenes(gene, reference_dict, coord_dict,vcf_dict):
    '''output: {strain:seq} dict for a geoi'''

    strain_seqs={}

    #get window start and stop for gene
    chrom=coord_dict[gene][1]
    window_start=coord_dict[gene][4]
    window_stop=coord_dict[gene][5]
    
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
    outfile=gene+'_amended.fa'
    n=80
    
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


#build chrom dict
print('parsing selscan outfiles...')
chrom_dict={}
for h12file in h12files:
    #get pos and vals
    pos_dict=ParseH12(h12file)
    #get chromosome names
    prefix=h12file.split('.')[0]
    chrom=prefix.split('chromosome')[1]
    if len(chrom)==1:
        chrom='chr0'+chrom
    else:
        chrom='chr'+chrom
    chrom_dict[chrom]=pos_dict

#get start and stop (and write outfile) of genes and their g1 windows
print('determining gene windows...')
coord_dict=getStartStop(ann_dict, chrom_dict) #spot checked for ESP1 (YGR098C), value corresponds to h12_h2h1 outfile

#get gene vars

print('loading variants in goi from vcf...')
vcf_dict=parseVCF(vcf_file, coord_dict)

#amend gene
print('amending genes...')
for gene in goi:
    amended_seqs, ref_seq=amendGenes(gene, reference_dict, coord_dict, vcf_dict)
    writeFasta(gene, amended_seqs, ref_seq)
    
