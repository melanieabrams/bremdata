import regex
import numpy as np
import sys
import subprocess as sp
import copy
import datetime
import time
import random

# PARAMETERS #    #FOR DIAGNOSING MERGE OF POOLED_READ FILES Version 3-19-2020 MBA, modified from cweiss's map_and_pool_BLAT.py and scorad rbseq
mapping_genome = "/usr2/people/carlyweiss/Amended_Genomes/Amended_Genomes/Concatenated/D1373_Z1.fa"  # location of the .udb database
min_ID = str(95) # minimum identity cutoff for mapping reads.  Must be a decimel proportion (0-1)
min_id_pooling = .95 # minimum identity for reads at the pooling stage.  Can be different then at mapping.  Must be (0-1)
gff = '/usr2/people/carlyweiss/SparScerinfo/YS2+CBS432+plasmid_clean' ## GFF to use for annotations.  
min_read_cutoff = 0  ## this is only for reporting summary statistics
min_genome_read_length = 40
min_len_align = 40
single_barcode = False #set to True if analyzing a test dataset with only one barcode
maskOffByOne = True # when True, this treats low abundance barcodes that are off-by-one as a sequencing error and removes them from analysis

min_id_pooling*=100

# HELP #l

if len(sys.argv) == 1:
    print("USAGE: python map_PB_reads.py out_directory inFile1._pooled_reads inFile2._pooled_reads")
    exit()
    
# INPUT # 
out_dir = sys.argv[1]
read_files = sys.argv[2:]  




# BEGIN FUNCTIONS #

def parse_pooled_reads(pooled_read_file): #function for beginning analysis *AFTER* pooled read outfiles already created for individual libraries
    split_loc_dict = {}  # will hold hierarchical data on each insertion site.  keys for nested dictionaries are scaffold, strand, position, barcode, value is # reads mapping there.

    f=open(pooled_read_file)
    next(f)
    for line in f:
        row_data=line.split('\t')
        barcode=row_data[0]
        chrom=row_data[1]
        strand=row_data[2]
        pos=int(row_data[3])
        annotation=row_data[4]
        number_of_reads=int(row_data[5])

        # add data to split_loc_dict #

        if chrom not in split_loc_dict:
            split_loc_dict[chrom] = {'+' : {}, '-' : {}}

        if pos not in split_loc_dict[chrom][strand]:
                split_loc_dict[chrom][strand][pos] = {barcode: number_of_reads}
                
        if barcode not in split_loc_dict[chrom][strand][pos]:
            split_loc_dict[chrom][strand][pos][barcode] = number_of_reads


    f.close()
    return split_loc_dict



def remove_multibarring (split_loc_dict, outfileprefix="fastq_filename_here"):
    '''removes barcodes that map to multiple independent insertions
    since those aren't useful for downstream analysis'''
    
    remaining_barcodes = [ ]
    duplicate_barcodes = {}
    
    #goes through dict and gets barcodes, finds duplicates
    for chrom in split_loc_dict:
        print("...discovering barcodes that map to multiple locations from "+chrom+"...")
        for strand in split_loc_dict[chrom]:
            for pos in split_loc_dict[chrom][strand]:
                for bc in split_loc_dict[chrom][strand][pos]:
                    if bc in remaining_barcodes:
##                        print(bc)
                        if bc in duplicate_barcodes:
                            duplicate_barcodes[bc]+=1
                        else:
                            duplicate_barcodes[bc]=2

                    else:
                        remaining_barcodes.append(bc)
                    
    unique_split_loc_dict = copy.deepcopy(split_loc_dict)
    
    #goes through dict and removes duplicate barcodes from further analysis
    for chrom in split_loc_dict:
        print("...discovering barcodes that map to multiple locations from "+chrom+"...")
        for strand in split_loc_dict[chrom]:
            for pos in split_loc_dict[chrom][strand]:
                for bc in split_loc_dict[chrom][strand][pos]:
                    if bc in duplicate_barcodes:
                        del unique_split_loc_dict[chrom][strand][pos][bc]
                        

    #writes separate outfile with duplicate barcodes
    wf =open(out_dir+outfileprefix+"_mulitmapping_barcodes", 'w')
    wf.writelines("barcode\t#of mapped locations\n")

    if len(duplicate_barcodes.keys()) > 0:
        for dupbar in duplicate_barcodes.keys():
            wf.writelines(dupbar+"\t"+str(duplicate_barcodes[dupbar])+"\n")
    wf.close()

    #writes number nonunique to mapping stats file
    number_nonunique = len(duplicate_barcodes)
    wf = open(out_dir+outfileprefix+"_mapping_stats", 'a')
    wf.writelines("\n\nRemoving non-unique barcodes: \n")
    wf.writelines("removed "+str(number_nonunique)+" non-unique barcode insertions\n")
    wf.close()


    print("removed "+str(number_nonunique)+" non-unique barcode insertions")
    return unique_split_loc_dict

def OffByOneList(seq):
    '''based on rbdnaseq pipeline, this
    generates a list of sequences off by one from the input
    use this in maskOffByOne to mask barcodes that are off by one due to sequencing error'''
    if seq[0] in ("A","T","G","C"):
        char_set = ("A","T","G","C")
    elif seq[0] in ("a","t","g","c"):
        char_set = ("a","t","g","c")
    else:
        return False

    variants = []
    for chari in list(range(len(seq))):
        if chari == 0:
            preseq = ""
        else:
            preseq = seq[0:chari]
        if chari == len(seq)-1:
            postseq = ""
        else:
            postseq = seq[chari+1:]
        for char in char_set:
            if seq[chari] != char:
                variants.append(preseq+char+postseq)
    return(variants)

def MaskOffByOne(split_loc_dict, outfileprefix="fastqfilename"):
    '''modified from rbdnaseq pipeline'''

    masked_split_loc_dict = copy.deepcopy(split_loc_dict)
    offByOneList = [0]
    totals = {}

    for chrom in split_loc_dict: #build a dictionary of barcode-to-total-read
        for strand in split_loc_dict[chrom]:
            for pos in split_loc_dict[chrom][strand]:
                for barcode in split_loc_dict[chrom][strand][pos]:
                    total = split_loc_dict[chrom][strand][pos][barcode]
                    if barcode in totals:
                        totals[barcode]+=total
                    else:
                        totals[barcode]=total
                    
    offByOneList = [ ]
    for chrom in split_loc_dict: #go through and remove barcodes where an off-by-one variant is 100 times more common
        print("...masking off-by-one barcodes in chromosome "+chrom+"...")
        for strand in split_loc_dict[chrom]:
            for pos in split_loc_dict[chrom][strand]:
                for barcode in split_loc_dict[chrom][strand][pos]:
                    variants = OffByOneList(barcode)
    ##                if barcode != 'TAAGCAACCTCGGCGCATAG':
    ##                    print(barcode)
    ##                    print(variants)
    ##                    print('TAAGCAACCTCGGCGCATAG' in variants)
    ##                    print(totals[barcode],totals['TAAGCAACCTCGGCGCATAG'])
                    offByOne = False
                    
                    try:
                        for variantBarcode in variants:
                            if (not offByOne) & (variantBarcode in totals):
                                if (totals[variantBarcode] > totals[barcode]*100):
                                    offByOneList.append(barcode)
                                    del masked_split_loc_dict[chrom][strand][pos][barcode]
                                    offByOne = True
                    except TypeError: #for some reason, one or two barcodes would throw a TypeError in one file if run in python2 (should be py3, but this guards against program aborting if run in old python). 
                        None

                
    print("masked "+str(len(offByOneList))+" off-by-one barcodes")
    wf = open(out_dir+outfileprefix+"_mapping_stats", 'a')
    wf.writelines("\n\nMasking off-by-one barcodes: \n")
    wf.writelines("masked "+str(len(offByOneList))+" off-by-one barcodes \n")
    wf.close()

    return masked_split_loc_dict 
    

def clean_barcodes(split_loc_dict, single_barcode=False, maskOffByOne=True, outfileprefix="read_file_name_here"):
    if single_barcode == False:
        split_loc_dict = remove_multibarring(split_loc_dict,outfileprefix=read_file)
       # print("...done removing non-unique barcodes...")
    if maskOffByOne == True:
        masked_unique_split_loc_dict = MaskOffByOne(split_loc_dict, outfileprefix=read_file)
        #print("...done masking off-by-one barcodes...")

    return split_loc_dict


def Annotate_insetions(split_loc_dict,read_file_name):

    out_filename = out_dir+read_file+"_pooled_reads"  # the final output, this will hold the pooled insertion table                                                                      
    wf = open(out_filename,'w')

    wf.writelines("ID\tscaffold\tstrand\tlocation\tannotation\tn\trel_loc\trel_prop\tgene_length\n")


    sc_insertions = 0
    sp_insertions = 0

    sc_genic_insertions = 0
    sp_genic_insertions = 0

    tot_insertions = 0
    tot_min_insertions = 0.0
    plasmid_reads = 0 # reads mapping to the plasmid
    Rtn_reads = 0 # reads mapping to the tn right border
    Ltn_reads = 0 # reads mapping to the tn left border

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


# Search through the dictionary for insertions that fall within genes                                                                                                                     
    for chrom in split_loc_dict:
        for strand in split_loc_dict[chrom]:
            for pos in split_loc_dict[chrom][strand]:
                for barcode in split_loc_dict[chrom][strand][pos]:
                    try:
                        if split_loc_dict[chrom][strand][pos][barcode] > min_read_cutoff:
                            tot_min_insertions+=1
                            if chrom[:2] == 'sc':
                                sc_insertions+=1
                            elif chrom[:2] == 'sp':
                                sp_insertions+=1
                        tot_insertions+=1
                    except TypeError:
                        print(split_loc_dict[chrom][strand][pos])
                        print(split_loc_dict[chrom][strand][pos][barcode])
                        print(min_read_cutoff)

                    # set defaults for noncoding                                                                                                                                                 

                    insertion_type = 'NC'
                    gene_length = -1
                    relative_insertion_site = -1
                    prop_gene_insertion_in = -1

                    for gff_chrom in gff_dict:
                        if gff_chrom != chrom:
                            continue
                        for gene in gff_dict[gff_chrom]:
                            if pos >= gff_dict[gff_chrom][gene][0] and pos <= gff_dict[gff_chrom][gene][1]:  # if the insertion falls within a gene                                              
                                insertion_type = gene
                                gene_length = gff_dict[gff_chrom][gene][1] - gff_dict[gff_chrom][gene][0]+1
                                if gff_dict[gff_chrom][gene][2] == '+':
                                    relative_insertion_site = pos - gff_dict[gff_chrom][gene][0]

                                else:
                                    relative_insertion_site = gff_dict[gff_chrom][gene][1] - pos

                                if gene[:2] == 'sp' and split_loc_dict[chrom][strand][pos][barcode] > min_read_cutoff:
                                    sp_genic_insertions+=1
                                elif gene[:2] == 'sc'  and split_loc_dict[chrom][strand][pos][barcode] > min_read_cutoff:
                                    sc_genic_insertions+=1

                                prop_gene_insertion_in = relative_insertion_site / float(gene_length)

                    wf.writelines(barcode+"\t"+chrom+"\t"+strand+"\t"+str(pos)+"\t"+insertion_type+"\t"+str(split_loc_dict[chrom][strand][pos][barcode])+"\t"+str(relative_insertion_site)+"\t"+str(prop_gene_insertion_in)+"\t"+str(gene_length)+"\n")

    wf.close()

    f = open(out_filename)
    for line in f:
        if line[:2] != 'pl':
            continue

        line = line.split("\t")
        if line[4] == 'TN_right_arm':
            Rtn_reads+=int(line[5]) # reads mapping to the tn right border        
        elif line[4] == 'TN_left_arm':

            Ltn_reads+=int(line[5]) # reads mapping to the tn right border        
        elif line[4] == 'NC':
            plasmid_reads+=int(line[5])

    f.close()

    
    tot_genic_insertions = float(sc_genic_insertions+sp_genic_insertions)

    wf = open(out_dir+read_file_name+"_mapping_stats", 'a')
    wf.writelines("\n\nAnnotating insertions: \n")
    wf.writelines("reads mapping to NC plasmid backbone: "+str(plasmid_reads)+"\n")
    wf.writelines("reads mapping to TN right border: "+str(Rtn_reads)+"\n")
    wf.writelines("reads mapping to TN left border: "+str(Ltn_reads)+"\n")

    wf.writelines("total insertions: "+str(tot_insertions)+"\n")
    wf.writelines("total insertions with >"+str(min_read_cutoff)+" reads: "+str(tot_min_insertions)+"\tscer: "+str(sc_insertions)+" ("+str(100*sc_insertions/tot_min_insertions)+"%)"+" spar: "+str(sp_insertions)+"("+str(100*sp_insertions/tot_min_insertions)+"%)\n")

    print("\ntotal insertions: "+str(tot_insertions))
    print("total genic insertions: "+str(tot_genic_insertions))

    if tot_genic_insertions > 0:
        
        wf.writelines("of these "+str(tot_insertions)+" insertions:\n")
        wf.writelines("total genic insertions: "+str(sp_genic_insertions + sc_genic_insertions)+" ("+str(100*(sp_genic_insertions + sc_genic_insertions)/tot_min_insertions)+"% of insertions)\n")
        wf.writelines("Scer genic insertions: "+str(sc_genic_insertions)+" ("+str(100*sc_genic_insertions/tot_genic_insertions)+"% of genic insertions)\n")
        wf.writelines("Spar genic insertions: "+str(sp_genic_insertions)+" ("+str(100*sp_genic_insertions/tot_genic_insertions)+"% of genic insertions)\n")

    else:
        wf.writelines("no genic insertions\n")


    wf.close()


def merge_all_tnseq(loc_dicts,merged_filename):
##    test_counter=0
##
##    test1=False
##    test2=False
##    test3=False
##
##    test_bc='TCCACAATGCACACTCCGTA'
##
##    print("test_bc in loc_dict[0]: "+str(test_bc in loc_dicts[0]))
##
          
    print("===merging_all_tnseq==")
    duplicate_reads = 0
    if len(loc_dicts)==1:
        return
    else:
        merged_dict = loc_dicts[0]
##        #test if in loc dict
##        test_in_0=False
##        for chrom in merged_dict:
##            for strand in merged_dict[chrom]:
##                for pos in merged_dict[chrom][strand]:
##                    for barcode in merged_dict[chrom][strand][pos]:
##                        if test_bc==barcode:
##                            test_in_0=True
##
##        print("test_bc in loc_dict[1]: "+str(test_in_0))
##
        #end test  if in loc_dict 
        duplicate_insert_dict = {}
##        if test_bc in duplicate_insert_dict:
##            test1=str(test_counter)+'_'+str(chrom)+':'+str(strand)+':'+str(pos)
##        test_in_1=False
        for loc_dict in loc_dicts[1:]:
##            test_counter+=1
            for chrom in loc_dict:
                if chrom not in merged_dict:
                    merged_dict[chrom] = {'+' : {}, '-' : {}} 
                for strand in loc_dict[chrom]:
                    for pos in loc_dict[chrom][strand]:
                        for barcode in loc_dict[chrom][strand][pos]:
##                            if test_bc==barcode:
##                                test_in_1=True
##                                print("hit barcode")
##                                print(str(pos in merged_dict[chrom][strand]))
##                                print(merged_dict[chrom][strand][pos])
                            total = loc_dict[chrom][strand][pos][barcode]
                            if pos not in merged_dict[chrom][strand]: # add an insert to the merged dict if not already one at that position
                                merged_dict[chrom][strand][pos]={}
                                merged_dict[chrom][strand][pos][barcode]=loc_dict[chrom][strand][pos][barcode]
                            elif barcode in merged_dict[chrom][strand][pos]: #if same barcode, add to total
##                                if barcode==test_bc:
##                                    test2=str(test_counter)+'_'+str(chrom)+':'+str(strand)+':'+str(pos)
##                                    print("test2: "+test2)
                                merged_dict[chrom][strand][pos][barcode]+=total
                                duplicate_reads+=1
                                if barcode in duplicate_insert_dict:
                                    duplicate_insert_dict[barcode][1]+=1
                                else:
                                    duplicate_insert_dict[barcode]=[chrom+strand+str(pos), 2]
                            else:
                                merged_dict[chrom][strand][pos][barcode] = loc_dict[chrom][strand][pos][barcode]
##                                if barcode==test_bc:
##                                    test3=str(test_counter)+'_'+str(chrom)+':'+str(strand)+':'+str(pos)
                                    
##        print("test_bc in loc_dict[1]: "+str(test_in_1))
##         
##        print("test barcode in duplicate_insert_dict?: "+str(test_bc in duplicate_insert_dict))                     
        wf=open(out_dir+merged_filename+"_duplicate_inserts", 'w')
        wf.writelines("combined "+str(duplicate_reads)+" duplicate inserts (same barcode at same position) from different fastq\n\n")
##        wf.writelines("print tests: 1="+str(test1)+", 2="+str(test2)+", 3="+str(test3)+" for bc "+test_bc+'\n')
        wf.writelines("barcode\tposition\tnumber of times found\n")
        for dupins in duplicate_insert_dict:
            wf.writelines(dupins+"\t"+duplicate_insert_dict[dupins][0]+"\t"+str(duplicate_insert_dict[dupins][1])+"\n")
            
        wf.close()
        
        
        merged_dict = clean_barcodes(merged_dict,single_barcode=single_barcode, maskOffByOne=maskOffByOne, outfileprefix=merged_filename)
        print("combined "+str(duplicate_reads)+" duplicate inserts (same barcode at same position) from different fastq")

        wf = open(out_dir+merged_filename+"_mapping_stats", 'a')
        wf.writelines('\n==='+merged_filename+'===\n')
        wf.writelines("BASH COMMAND\n"+command_line_string+ "\nRUN DATETIME\n" + now.strftime("%m/%d/%Y, %H:%M:%S")+"\n")
        wf.writelines("\nMerging libraries:\n")
        wf.writelines("combined "+str(duplicate_reads)+" duplicate inserts (same barcode at same position) from different fastq")
        wf.close()

        
        
        Annotate_insetions(merged_dict, merged_filename) # identify insertions in genes, and write the final pooled outfile for a given fastq
        print ("...done annotating")
    return
                            
                        


                
#### START PROGRAM ####

command_line_string = "" #get command line so can add to outfiles
for argv in sys.argv[:]:
    command_line_string+=(str(argv)+" ")

start_time = time.time() #Time before the operations start'
now = datetime.datetime.now()

loc_dicts = [] 
for read_file in read_files:
    print('...parsing'+read_file)
    loc_dict = parse_pooled_reads(read_file) #parse from already written file
    loc_dicts.append(loc_dict)

if len(loc_dicts)>1:
    print()
    merged_filename = "merged_tnseq_"+now.strftime("%Y-%m-%d")+'.'
    merge_all_tnseq(loc_dicts,merged_filename) #merge dictionaries and make final pooled outfile for all fastq

end_time = time.time() #Time after it finished

print("Merge took ", end_time-start_time, " seconds")
