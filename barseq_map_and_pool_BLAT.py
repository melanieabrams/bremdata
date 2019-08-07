import regex
import numpy as np
import sys
import subprocess as sp
import copy
import datetime

# PARAMETERS #    #Version 8-6-19 MBA, modified from cweiss's map_and_pool_BLAT.py and scorad rbseq

mapping_genome = "/usr2/people/carlyweiss/Amended_Genomes/Amended_Genomes/Concatenated/D1373_Z1.fa"  # location of the .udb database
min_ID = str(95) # minimum identity cutoff for mapping reads.  Must be a decimel proportion (0-1)
min_id_pooling = .95 # minimum identity for reads at the pooling stage.  Can be different then at mapping.  Must be (0-1)
gff = '/usr2/people/carlyweiss/SparScerinfo/YS2+CBS432+plasmid_clean' ## GFF to use for annotations.  
min_read_cutoff = 0  ## this is only for reporting summary statistics
min_genome_read_length = 40
min_len_align = 40
single_barcode = True #set to True if analyzing a test dataset with only one barcode
maskOffByOne = True # when True, this treats low abundance barcodes that are off-by-one as a sequencing error and removes them from analysis

min_id_pooling*=100

# HELP #

if len(sys.argv) == 1:
    print("USAGE: python map_PB_reads.py out_directory fastq_file1 fastq_file2")
    exit()
    
# INPUT # 
out_dir = sys.argv[1]
read_files = sys.argv[2:]  




# BEGIN FUNCTIONS # 

def Filter_Reads(fastq_file):
    '''modified from rh-seq, filter out reads that do not have a TN sequence or the universal primers for barcodes'''

    wf = open(out_dir+fastq_filename+'_parsed_reads','w')  # outfile for the trucated reads that will be mapped

    #organization: ....flanking_left-->barcode-->flanking_right-->...-->tn_patten-->genome_read-->rd2_sequence_primer
    tn_pattern = regex.compile('(?e)(ATGATTATCTTTCTAGGGTTAA){e<=2}')  # last 22 bps of the right arm of the transposon.  searches for this pattern allowing 2 mismatchs for sequencing errors
    rd2_pattern = regex.compile('(?e)(AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC){e<=2}')  # Universal primer 1.  searches for this pattern allowing 2 mismatchs for sequencing errors
    
    flanking_left_pattern = regex.compile('(?e)(GATGTCCACGAGGTCTCT){e<=2}')  # Universal primer 1.  searches for this pattern allowing 2 mismatchs for sequencing errors
    flanking_right_pattern = regex.compile('(?e)(CGTACGCTGCAGGTCGAC){e<=2}')  # Universal primer 2.  searches for this pattern allowing 2 mismatchs for sequencing errors
    # counts for summary stats #

    tot_reads = 0.0  
    reads_with_tn = 0
    line_count = 0
    total_Ns = 0
    total_bases = 0
    reads_too_short = 0
    too_high_Ns = 0
    max_length_genome_read = 0 
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

            ###make sure transposon with barcode primers hopped in to read of interest
            fl_match_data = flanking_left_pattern.search(read)  # searches for the tn pattern
            if fl_match_data == None:
                continue  #skips the read if it doesn't have universal primer 1

            fr_match_data = flanking_right_pattern.search(read)  # searches for the tn pattern
            if fr_match_data == None:
                continue  #skips the read if it doesn't have a tn pattern
            
            tn_match_data = tn_pattern.search(read)  # searches for the tn pattern
            if tn_match_data == None:
                continue  #skips the read if it doesn't have a tn pattern

            rd2_match_data = rd2_pattern.search(read)  # searches for the tn pattern
            if rd2_match_data == None:
                genome_read = read[tn_match_data.end()-4:] # gets the sequence of the read outside the transposon (still contains TTAA)
            else:
                genome_read =  read[tn_match_data.end()-4:rd2_match_data.start()] #also removes end primer if present

            barcode = read[fl_match_data.end():fr_match_data.start()]

            
            #calculate stats
            total_bases+=len(read)
            total_Ns+=read.count('N') # counts the numer of Ns per read
            reads_with_tn+=1  ## counts reads that have the tn sequence
            end_match = tn_match_data.end()  # finds the location in the read that the last tn base matches to
            
            longest_read = '' #update max genome length read...useful for diagnosis if all of them are too short
            if len(genome_read) > max_length_genome_read:
                max_length_genome_read = len(genome_read)
                longest_read = read
                longest_genome_read = genome_read
    

            #throw away reads that are too short or have too many Ns
            if len(genome_read) < min_genome_read_length: # if the length of the parsed read is less than min, throw it out.  50bp here is appropriate for 150bp reads.
                reads_too_short+=1
                continue
            
            if genome_read.count('N') / float(len(genome_read)) > .2:
                too_high_Ns+=1
                continue

            
            # writes a new fastq file that has the filtered read and only the regions TTAA--genome.  The TN sequence is removed, and the header is appended with the barcode
            wf.writelines(">"+barcode+"\n")
            wf.writelines(genome_read+"\n")
           # wf.writelines("+\n")
           # wf.writelines(qual[tn_match_data.end()-4:])
            
    wf.close()    

    prop_Ns = float(total_Ns) / total_bases # gets the proportion of Ns in the fastq file as an indicator for overall sequence quality
    prop_Ns = round(prop_Ns,2)

    wf = open(out_dir+fastq_filename+"_mapping_stats", 'w')  # the file that will contain summary stats of the mapping and additional analysis

    # write mapping statistics" 

    print("total_reads: ", str(tot_reads))
    wf.writelines("total_reads: "+str(tot_reads)+"\n")
    print("reads with tn: ", str(reads_with_tn), " ("+str(100*reads_with_tn/tot_reads)[:4]+"%)")
    wf.writelines("reads with TN: "+str(reads_with_tn)+" ("+str(100*reads_with_tn/tot_reads)[:4]+"%)\n")
    print("proportion bases in TN containing reads that are 'N': ", str(prop_Ns))
    wf.writelines("proportion bases in TN containing reads that are 'N': "+ str(prop_Ns)+"\n")
    print("TN contaning reads too short that were discarded: ", str(reads_too_short), " ("+str(100*reads_too_short/float(reads_with_tn))[:4]+"%)")
    wf.writelines("TN contaning reads too short that were discarded: "+str(reads_too_short)+" ("+str(100*reads_too_short/float(reads_with_tn))[:4]+"%)\n")
    print("max genome read length: " + str(max_length_genome_read))
##    print("longest_read:")
##    print(longest_read)
##    print(longest_genome_read)
##    print(len(longest_genome_read) < min_genome_read_length)
    wf.writelines("max genome read length: " + str(max_length_genome_read))
    print("TN containing reads with too many Ns: ", str(too_high_Ns), " ("+str(100*too_high_Ns/float(reads_with_tn))+"%)")
    wf.writelines("TN containing reads with too many Ns: "+str(too_high_Ns)+" ("+str(100*too_high_Ns/float(reads_with_tn))[:4]+"%)\n")
    print("barcodes containing reads with too many Ns: ", str(too_high_Ns), " ("+str(100*too_high_Ns/float(reads_with_tn))[:4]+"%)")
    wf.writelines("TN containing reads with too many Ns: "+str(too_high_Ns)+" ("+str(100*too_high_Ns/float(reads_with_tn))[:4]+"%)\n")

    wf.close()

    tot_parsed_reads = reads_with_tn - reads_too_short - too_high_Ns

    return tot_parsed_reads


def Map_Reads():
    ''' from rh-seq, uses blat to map the reads'''

#    cmd = ["/opt/bin/usearch8.1.1756_i86linux64", "-threads", nthreads, "-usearch_local",  out_dir+fastq_filename+"_parsed_reads", "-db", db, "-strand", "both", "-id", min_id, "-alnout", out_dir+fastq_filename+"_aligned_reads", "-maxhits", "2", "-notmatched", out_dir+fastq_filename+"_unmapped_reads", "-userout", out_dir+fastq_filename+"_mapped_reads", "-userfields", "query+target+id+alnlen+qstrand+tstrand+tlot+thit", "-mincols", "50", "-query_cov", ".95"]  # the command used to call usearch

    cmd = ["blat", mapping_genome, out_dir+fastq_filename+"_parsed_reads", out_dir+fastq_filename+"_mapped_reads", "-minIdentity=95", "-tileSize=12", "-out=blast8" ]
    
    sp.call(cmd)

def Filter_for_ID_and_multimapping_and_pool(num_parsed_reads):
    ''' modified from rh-seq
    filters our mapped reads that map below the ID threshold, those that map to multiple locations, then pools reads into insertion sites
    also filters out reads with barcodes that show up more than once'''

    mapped_file = out_dir+fastq_filename+"_mapped_reads" # where the mapped reads are

    first = 'y'

    f = open(mapped_file)


    read_dict = {}
    total_mapped_reads = 0 ###  CHANGE THIS BACK TO 0
    reads_above_identity_threshold = 0
    mulitmapped_reads = 0

    all_read_list = set()

    for line in f:

        
        line = line.strip().split("\t")

        
        read = line[0]
       
        if read not in all_read_list:  ## counts all reads that map at least once
            all_read_list.add(read)
            total_mapped_reads+=1

        ID = float(line[2])
        if ID < min_id_pooling:  ## filters out mapped reads below the minimum identity for pooling                                                                                                
            continue

        len_align = float(line[3])
        start_align = float(line[6]) 

        if len_align < min_len_align or start_align > 3:
            continue

        reads_above_identity_threshold+=1  # this counts reads that are above the identity threshold and length of aligment threshold and in which the alignment starts within 3bps of the read TTAA
        
        start = line[8]
        end = line[9]

        strand = '+'
        if start < end:
            strand = '-'
        

        insertion_loc = start
        
        loc = line[1]+"__"+strand+"__"+insertion_loc  ## This is the identifier for the insertion. scaffold+strand+position


        # searches for reads that map to more than 1 locatiion # 
        if read in read_dict:  
            read_dict[read].append([loc, ID])  ## appends the read location and %ID of mapping in case of multiple mappings of the same read
        #    mulitmapped_reads+=2  # as long as the usearch multimapping paramter is set to a maximim of 2 reads, this will work.  
        else:
            read_dict[read] = [[loc,ID]]


    # filters our reads that map to more than one location # 


    print("...done making the read_dict...")
    
    barcode_lengths = {}
    for barcode in read_dict.keys():
        if len(barcode) in barcode_lengths.keys():
            barcode_lengths[len(barcode)]+=1
        else:
            barcode_lengths[len(barcode)]=1
        
    loc_dict = {}  # 
    for barcode in read_dict.keys():
        for map_loc in range(len(read_dict[barcode])):
            loc_name=read_dict[barcode][map_loc][0]
            if loc_name in loc_dict:
                loc_dict[loc_name][1]+=1
            else:
                loc_dict[loc_name]=[barcode,1]

    print("total_mapped barcodes: ", str(total_mapped_reads), " ("+str(float(100*total_mapped_reads/num_parsed_reads))[:4]+"% of parsed reads)")
    print("mapped reads passing identity cutoff: "+str(reads_above_identity_threshold), "("+str(100*reads_above_identity_threshold/float(total_mapped_reads))[:4]+"% of mapped barcodes)")
    print("remaining reads mapping to one location: "+str(len(read_dict)), "("+str(100*len(read_dict)/float(reads_above_identity_threshold))[:4]+"%)")

    wf = open(out_dir+fastq_filename+"_mapping_stats", 'a')

    wf.writelines("total mapped reads: "+str(total_mapped_reads)+" ("+str(float(100*total_mapped_reads/num_parsed_reads))[:4]+"% of parsed reads)\n")
    wf.writelines("mapped reads passing identity cutoff: "+str(reads_above_identity_threshold)+" "+str(100*reads_above_identity_threshold/float(total_mapped_reads))[:4]+"% of mapped barcodes\n")
    wf.writelines("remaining barcodes mapping to one location: "+str(len(read_dict))+" "+str(100*len(read_dict)/float(reads_above_identity_threshold))[:4]+"%\n")

    wf.close()

    return loc_dict, len(read_dict)


def Combine_near_mappings(loc_dict, reads_remaining):
    '''modified from rh-seq, combines insertion sites that are within 3 bases of each other.  reads are assigned to the site with the initial max number of reads'''

    split_loc_dict = {}  # will hold hierarchical data on each insertion site.  keys for nested dictionaries are scaffold, strand, position and value is # reads mapping there. 

    for full_location in loc_dict:  ## loc dict holds the identifier for an inserion site as key and reads mapping to that site as a value
        barcode = loc_dict[full_location][0]
        number_of_reads = loc_dict[full_location][1]
        chrom = full_location.split("__")[0]
        strand = full_location.split("__")[1]
        pos = int(full_location.split("__")[2])

        # initialize the dictionary #

        if chrom not in split_loc_dict:
            split_loc_dict[chrom] = {'+' : {}, '-' : {}}

        if pos not in split_loc_dict[chrom][strand]:
            split_loc_dict[chrom][strand][pos] = [number_of_reads,barcode]

    reads_moved = 0
    # sorts the insertion positions, and combines reads forward, then reverses the sorting and combines forward again.  #

    for chrom in split_loc_dict:
        for strand in split_loc_dict[chrom]:

            sorted_positions = sorted(split_loc_dict[chrom][strand])
            first ='y'
            for pos in sorted_positions:
                if first == 'y':
                    first = 'n'
                    last = pos
                    continue

                if int(pos) - int(last) < 4:

                    if split_loc_dict[chrom][strand][pos][0] >= split_loc_dict[chrom][strand][last][0]:
                        split_loc_dict[chrom][strand][pos][0]+=split_loc_dict[chrom][strand][last][0]
                        reads_moved+=split_loc_dict[chrom][strand][last][0]
                        del split_loc_dict[chrom][strand][last]
                last = pos
                
            sorted_positions = sorted(split_loc_dict[chrom][strand])
            sorted_positions.reverse()

            first ='y'
            for pos in sorted_positions:
                if first == 'y':
                    first = 'n'
                    last = pos
                    continue

                if abs(int(pos) - int(last)) < 4:
                    if split_loc_dict[chrom][strand][pos][0] >= split_loc_dict[chrom][strand][last][0]:
                        split_loc_dict[chrom][strand][pos][0]+=split_loc_dict[chrom][strand][last][0]
                        reads_moved+=split_loc_dict[chrom][strand][last][0]
                        del split_loc_dict[chrom][strand][last]

                last = pos
    print("remaining reads moved to higher peak: ", str(reads_moved), "("+str(100*float(reads_moved)/reads_remaining)[:4]+"% of mapped barcodes)")
    

    wf = open(out_dir+fastq_filename+"_mapping_stats", 'a')
    wf.writelines("remaining reads moved to higher peak: "+str(reads_moved)+" "+str(100*float(reads_moved)/reads_remaining)[:4]+"% of mapped barcodes\n")
    wf.close()

    return split_loc_dict

def remove_multibarring (split_loc_dict, mapped_reads):
    '''removes barcodes that map to multiple independent insertions
    since those aren't useful for downstream analysis'''
    remaining_barcodes = [ ]
    duplicate_barcodes = [ ]
    #goes through dict and gets barcodes, finds duplicates
    for chrom in split_loc_dict:
        for strand in split_loc_dict[chrom]:
            for pos in split_loc_dict[chrom][strand]:
                bc = split_loc_dict[chrom][strand][pos][1]
                if bc in remaining_barcodes:
                    duplicate_barcodes.append(bc)
                else:
                    remaining_barcodes.append(bc)
                    
    unique_split_loc_dict = copy.deepcopy(split_loc_dict)                 
    #goes through dict and removes duplicate barcodes
    for chrom in split_loc_dict:
            for strand in split_loc_dict[chrom]:
                for pos in split_loc_dict[chrom][strand]:
                    bc = split_loc_dict[chrom][strand][pos][1]
                    if bc in duplicate_barcodes:
                        del unique_split_loc_dict[chrom][strand][pos]
    number_nonunique = len(duplicate_barcodes)
    print("removed "+str(number_nonunique)+" non-unique barcode insertions")
    wf = open(out_dir+fastq_filename+"_mapping_stats", 'a')
    wf.writelines("removed "+str(number_nonunique)+" non-unique barcode insertions")
    return unique_split_loc_dict

def OffByOneList(seq):
    '''from rbdnaseq pipeline
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

def MaskOffByOne(split_loc_dict, mapped_reads):
    '''modified from rbdnaseq pipeline'''

    masked_split_loc_dict = copy.deepcopy(split_loc_dict)
    offByOneList = [0]
    totals = {}

    for chrom in split_loc_dict: #build a dictionary of barcode-to-total-read
        for strand in split_loc_dict[chrom]:
            for pos in split_loc_dict[chrom][strand]:
                barcode = split_loc_dict[chrom][strand][pos][1]
                total = split_loc_dict[chrom][strand][pos][0]
                if barcode in totals:
                    totals[barcode]+=total
                else:
                    totals[barcode]=total
                    
    offByOneList = [ ]
    for chrom in split_loc_dict: #go through and remove barcodes where an off-by-one variant is 100 times more common
        for strand in split_loc_dict[chrom]:
            for pos in split_loc_dict[chrom][strand]:
                barcode = split_loc_dict[chrom][strand][pos][1]
                variants = OffByOneList(barcode)
##                if barcode != 'TAAGCAACCTCGGCGCATAG':
##                    print(barcode)
##                    print(variants)
##                    print('TAAGCAACCTCGGCGCATAG' in variants)
##                    print(totals[barcode],totals['TAAGCAACCTCGGCGCATAG'])
                offByOne = False
                for variantBarcode in variants:
                    if (not offByOne) & (variantBarcode in totals):
                        if (totals[variantBarcode] > totals[barcode]*100):
                            offByOneList.append(barcode)
                            del masked_split_loc_dict[chrom][strand][pos]
                            offByOne = True

                
    print("masked "+str(len(offByOneList))+" off-by-one barcodes")
    return masked_split_loc_dict 
    

def clean_barcodes(split_loc_dict, mapped_reads, single_barcode=False, maskOffByOne=True):
    if single_barcode == False:
        split_loc_dict = remove_multibarring(split_loc_dict, mapped_reads)
        print("...done removing non-unique barcodes...")
    if maskOffByOne == True:
        masked_unique_split_loc_dict = MaskOffByOne(split_loc_dict, mapped_reads)
        print("...done masking off-by-one barcodes...")
    return split_loc_dict


def Annotate_insetions(split_loc_dict, mapped_reads,fastq_filename):

    out_filename = out_dir+fastq_filename+"_pooled_reads"  # the final output, this will hold the pooled insertion table                                                                      
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
                if split_loc_dict[chrom][strand][pos][0] > min_read_cutoff:
                    tot_min_insertions+=1
                    if chrom[:2] == 'sc':
                        sc_insertions+=1
                    elif chrom[:2] == 'sp':
                        sp_insertions+=1
                tot_insertions+=1

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

                            if gene[:2] == 'sp' and split_loc_dict[chrom][strand][pos][0] > min_read_cutoff:
                                sp_genic_insertions+=1
                            elif gene[:2] == 'sc'  and split_loc_dict[chrom][strand][pos][0] > min_read_cutoff:
                                sc_genic_insertions+=1

                            prop_gene_insertion_in = relative_insertion_site / float(gene_length)

                wf.writelines(str(split_loc_dict[chrom][strand][pos][1])+"\t"+chrom+"\t"+strand+"\t"+str(pos)+"\t"+insertion_type+"\t"+str(split_loc_dict[chrom][strand][pos][0])+"\t"+str(relative_insertion_site)+"\t"+str(prop_gene_insertion_in)+"\t"+str(gene_length)+"\n")

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

    wf = open(out_dir+fastq_filename+"_mapping_stats", 'a')

    wf.writelines("reads mapping to NC plasmid backbone: "+str(plasmid_reads)+" ("+str(100*plasmid_reads/mapped_reads)+"%) of mapped reads\n")
    wf.writelines("reads mapping to TN right border: "+str(Rtn_reads)+" ("+str(100*Rtn_reads/mapped_reads)+"%) of mapped reads\n")
    wf.writelines("reads mapping to TN left border: "+str(Ltn_reads)+" ("+str(100*Ltn_reads/mapped_reads)+"%) of mapped reads\n")

    wf.writelines("total insertions: "+str(tot_insertions)+"\n")
    wf.writelines("total insertions with >"+str(min_read_cutoff)+" reads: "+str(tot_min_insertions)+"\tscer: "+str(sc_insertions)+" ("+str(100*sc_insertions/tot_min_insertions)+"%)"+" spar: "+str(sp_insertions)+"("+str(100*sp_insertions/tot_min_insertions)+"%)\n")


    if tot_genic_insertions > 0:
        


        wf.writelines("OF  THESE:\n")
        wf.writelines("total genic insertions: "+str(sp_genic_insertions + sc_genic_insertions)+" ("+str(100*(sp_genic_insertions + sc_genic_insertions)/tot_min_insertions)+"% of insertions)\n")
        wf.writelines("Scer genic insertions: "+str(sc_genic_insertions)+" ("+str(100*sc_genic_insertions/tot_genic_insertions)+"% of genic insertions)\n")
        wf.writelines("Spar genic insertions: "+str(sp_genic_insertions)+" ("+str(100*sp_genic_insertions/tot_genic_insertions)+"% of genic insertions)\n")

    else:
        wf.writelines("no genic insertions\n")


    wf.close()


def merge_all_tnseq(loc_dicts,merged_filename):
    print("merging_all_tnseq")
    duplicate_reads = 0
    if len(loc_dicts)==1:
        return
    else:
        merged_dict = loc_dicts[0]
        for loc_dict in loc_dicts[1:]:
            for chrom in loc_dict:
                if chrom not in merged_dict:
                    merged_dict[chrom] = {'+' : {}, '-' : {}}
                for strand in loc_dict[chrom]:
                    for pos in loc_dict[chrom][strand]:
                        barcode = loc_dict[chrom][strand][pos][1]
                        total = loc_dict[chrom][strand][pos][0]
                        if pos not in merged_dict[chrom][strand]: # add an insert to the merged dict if not already one at that position
                            merged_dict[chrom][strand][pos]=loc_dict[chrom][strand][pos]
                        elif barcode == merged_dict[chrom][strand][pos][1]: #if same barcode, add to total
                            merged_dict[chrom][strand][pos][0]+=total
                            duplicate_reads+=1
                        else:
                            merged_dict[chrom][strand][pos].append(loc_dict[chrom][strand][pos])
        print("combined "+str(duplicate_reads)+" duplicate inserts from different fastq")
        merged_dict = clean_barcodes(merged_dict, reads_remaining,single_barcode=single_barcode, maskOffByOne=maskOffByOne)
        Annotate_insetions(merged_dict, reads_remaining, merged_filename) # identify insertions in genes, and write the final pooled outfile for a given fastq
        print ("...done annotating")
    return
                            
                        


                
#### START PROGRAM ####

loc_dicts = [] 

for read_file in read_files:
    fastq_filename = read_file.split("/")[-1] #get file name 
    print('\n'+fastq_filename+'\n') 
    num_parsed_reads = Filter_Reads(read_file)  ## filters out reads that don't have tn sequence, writes the genomic portion of the remaining reads to a new file
    Map_Reads()  ## maps reads
    print("...done mapping...")
    loc_dict, reads_remaining = Filter_for_ID_and_multimapping_and_pool(num_parsed_reads)  # filters out reads below the identity threshold, that map to multiple locations and then pools insertions 
    print("...done filtering...")
    loc_dict = Combine_near_mappings(loc_dict, reads_remaining) # combine insertions within 3bp of each other
    print("...done combine near mappings...")
    loc_dict = clean_barcodes(loc_dict, reads_remaining,single_barcode=single_barcode, maskOffByOne=maskOffByOne)
    Annotate_insetions(loc_dict, reads_remaining, fastq_filename) # identify insertions in genes, and write the final pooled outfile for a given fastq
    print ("...done annotating")
    loc_dicts.append(loc_dict)

if len(loc_dicts)>1:
    now = datetime.datetime.now()
    print()
    merged_filename = "merged_tnseq_"+now.strftime("%Y-%m-%d")+'.'
    merge_all_tnseq(loc_dicts,merged_filename) #merge dictionaries and make final pooled outfile for all fastq
