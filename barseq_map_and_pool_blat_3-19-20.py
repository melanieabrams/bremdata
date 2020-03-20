import regex
import numpy as np
import sys
import subprocess as sp
import copy
import datetime
import time
import random

# PARAMETERS #    #Version 2-19-2020 MBA, modified from cweiss's map_and_pool_BLAT.py and scorad rbseq

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


    read_dict = {}  #format of read dict will be {barcode:[chrom_strand_pos, ID]}
    total_mapped_barcodes = 0 
    reads_above_identity_threshold = 0
    mulitmapped_reads = 0

    all_barcode_list = set()
    number_repeat_barcodes = 0

    for line in f:
        line = line.strip().split("\t")
        bc = line[0]       
        if bc not in all_barcode_list:  ## this counts the number of barcodes that mapped to at least one location in the genome
            all_barcode_list.add(bc)
            total_mapped_barcodes+=1
        else:
            number_repeat_barcodes+=1

      
        ID = float(line[2])
  
        if ID < min_id_pooling:  ## filters out mapped reads below the minimum mapping identity                                                                                               
            continue

        len_align = float(line[3])
        start_align = float(line[6]) 
        if len_align < min_len_align or start_align > 3: # filters out mapped reads that are below the length of alignment threshold or where the alignment starts more than 3bps after the TTAA
            continue

        reads_above_identity_threshold+=1  # this counts reads that are above the identity threshold and length of aligment threshold and in which the alignment starts within 3bps of the read TTAA
        
        start = line[8]
        end = line[9]

        strand = '+'
        if start < end:
            strand = '-'
        
        insertion_loc = start
        
        loc = line[1]+"__"+strand+"__"+insertion_loc  ## This is an insertion's scaffold+strand+position.  It used to be the identifier for an inset before barcodes were used in RH-Seq.


        # searches for reads that map to more than 1 locatiion # 
        if bc in read_dict:  
            read_dict[bc].append([loc, ID])  ## appends the read location and %ID of mapping in case of multiple mappings of the same read
        #    mulitmapped_reads+=2  # as long as the usearch multimapping paramter is set to a maximim of 2 reads, this will work.  
        else:
            read_dict[bc] = [[loc,ID]]


    print("...done making the read_dict...")
    
    barcode_lengths = {} #collects stats on how long the barcodes in the library were
    for barcode in read_dict.keys():
        if len(barcode) in barcode_lengths.keys():
            barcode_lengths[len(barcode)]+=1
        else:
            barcode_lengths[len(barcode)]=1
        
    loc_dict = {}  #format of loc_dict is {position:{barcode: occurences}}
    number_occurences=0
    number_of_positions_with_multiple_barcodes=0
    for barcode in read_dict.keys():
        for map_loc in range(len(read_dict[barcode])):
            loc_name=read_dict[barcode][map_loc][0]
            if loc_name in loc_dict:
                if barcode in loc_dict[loc_name].keys():
                    loc_dict[loc_name][barcode]+=1
                else:
                    loc_dict[loc_name][barcode]=1
                    number_of_positions_with_multiple_barcodes=0
            else:
                loc_dict[loc_name]={barcode:1}



    print("out of the total reads parsed from the blat mapping outfile, there were "+str(total_mapped_barcodes+number_repeat_barcodes)+" reads with barcodes")
    print("--total number of different barcodes that were mapped: " + str(total_mapped_barcodes)+" ("+ str(100*total_mapped_barcodes/(total_mapped_barcodes+number_repeat_barcodes))[:4]+"% of parsed reads with barcodes )")
    print("--number of times the script finds a read containing a barcode it had already encountered: "+str(number_repeat_barcodes)+" ("+ str(100*number_repeat_barcodes/(total_mapped_barcodes+number_repeat_barcodes))[:4]+"% of parsed reads with barcodes)\n")
    print("...now, doing some filtering...")
    print("reads where the genomic portion passes mapping identity cutoff: "+str(reads_above_identity_threshold)+" ("+ str(100*reads_above_identity_threshold/(total_mapped_barcodes+number_repeat_barcodes))[:4]+"% of parsed reads with barcodes )")
    print("--remaining barcodes mapping to a genomic location: "+str(len(read_dict))+" ("+ str(100*len(read_dict)/total_mapped_barcodes)[:4]+"% of barcodes encountered)")
    print("genomic locations where more than one barcode inserted: "+str(number_of_positions_with_multiple_barcodes)+"\n")

    wf = open(out_dir+fastq_filename+"_mapping_stats", 'a') 

    wf.writelines("\nout of the total reads parsed from the blat mapping outfile, there were "+str(total_mapped_barcodes+number_repeat_barcodes)+" reads with barcodes\n")
    wf.writelines("total number of different barcodes that were mapped: " + str(total_mapped_barcodes)+" ("+ str(100*total_mapped_barcodes/(total_mapped_barcodes+number_repeat_barcodes))[:4]+"% of parsed reads with barcodes)\n")
    wf.writelines("number of times the script encounters a barcode it had already seen: "+str(number_repeat_barcodes)+" ("+ str(100*number_repeat_barcodes/(total_mapped_barcodes+number_repeat_barcodes))[:4]+"% of parsed reads with barcodes)\n")
    wf.writelines("reads where the genomic portion passed mapping identity cutoff: "+str(reads_above_identity_threshold)+" ("+ str(100*reads_above_identity_threshold/(total_mapped_barcodes+number_repeat_barcodes))[:4]+"% of parsed reads with barcodes)\n")
    wf.writelines("remaining barcodes mapping to a location: "+str(len(read_dict))+str(len(read_dict))+" ("+ str(100*len(read_dict)/total_mapped_barcodes)[:4]+"% of barcodes encountered)\n")
    wf.writelines("genomic locations where more than one barcode inserted: "+str(number_of_positions_with_multiple_barcodes))


    wf.close()


    return loc_dict, len(read_dict)



def Combine_near_mappings(loc_dict, reads_remaining):
    '''modified from rh-seq, combines insertion sites that are within 3 bases of each other.  reads are assigned to the site with the initial max number of reads'''    

    split_loc_dict = {}  # will hold hierarchical data on each insertion site.  keys for nested dictionaries are scaffold, strand, position, barcode, value is # reads mapping there.
    for full_location in loc_dict:  ## loc dict holds the identifier for an inserion site as key and reads mapping to that site as a value
        for barcode in loc_dict[full_location]:
            number_of_reads = loc_dict[full_location][barcode]
            chrom = full_location.split("__")[0]
            strand = full_location.split("__")[1]
            pos = int(full_location.split("__")[2])

            # initialize the dictionary #

            if chrom not in split_loc_dict:
                split_loc_dict[chrom] = {'+' : {}, '-' : {}}

            if pos not in split_loc_dict[chrom][strand]:
                split_loc_dict[chrom][strand][pos] = {barcode: number_of_reads}

            if barcode not in split_loc_dict[chrom][strand][pos]:
                split_loc_dict[chrom][strand][pos][barcode] = number_of_reads
            
    # sorts the insertion positions, and combines reads forward, then reverses the sorting and combines forward again.  #
    reads_moved = 0

    for chrom in split_loc_dict:
        for strand in split_loc_dict[chrom]:

            sorted_positions = sorted(split_loc_dict[chrom][strand]) #for each chromosome, and for each strand, sorts position
            first ='y'
            for pos in sorted_positions: #then for each of the positions on a given strand on a given chromosome
                if first == 'y': #skips the first because it can't combine with another yet, till combines backwards
                    first = 'n'
                    last = pos
                    continue
            

                if int(pos) - int(last) < 4:
                    delete_from_last=[]
                    for barcode in split_loc_dict[chrom][strand][last]:
                        if barcode in split_loc_dict[chrom][strand][pos]:
                            if split_loc_dict[chrom][strand][pos][barcode]>= split_loc_dict[chrom][strand][last][barcode]:
                                split_loc_dict[chrom][strand][pos][barcode]+=split_loc_dict[chrom][strand][last][barcode]
                                reads_moved+=split_loc_dict[chrom][strand][last][barcode]
                                delete_from_last.append(barcode) #remove barcodes whose reads moved to a higher peak from 'last'
                    for barcode in delete_from_last:
                        del split_loc_dict[chrom][strand][last][barcode]
                       
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
                    delete_from_last=[]
                    for barcode in split_loc_dict[chrom][strand][last]:
                        if barcode in split_loc_dict[chrom][strand][pos]:
                            if split_loc_dict[chrom][strand][pos][barcode] >= split_loc_dict[chrom][strand][last][barcode]:
                                split_loc_dict[chrom][strand][pos][barcode]+=split_loc_dict[chrom][strand][last][barcode]
                                reads_moved+=split_loc_dict[chrom][strand][last][barcode]
                                delete_from_last.append(barcode)
                    for barcode in delete_from_last:
                        del split_loc_dict[chrom][strand][last][barcode]

                last = pos
    print("remaining reads moved to higher peak when script combines near mappings: ", str(reads_moved))
    

    wf = open(out_dir+fastq_filename+"_mapping_stats", 'a')
    wf.writelines("\n\nCombining near mappings: \n")
    wf.writelines("remaining reads moved to higher peak: "+str(reads_moved)+"\n")
    wf.close()

    return split_loc_dict

def remove_multibarring (split_loc_dict, mapped_reads,outfileprefix="fastq_filename_here"):
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

def MaskOffByOne(split_loc_dict, mapped_reads, outfileprefix="fastqfilename"):
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
    

def clean_barcodes(split_loc_dict, mapped_reads, single_barcode=False, maskOffByOne=True, outfileprefix="fastq_filename_here"):
    if single_barcode == False:
        split_loc_dict = remove_multibarring(split_loc_dict, mapped_reads,outfileprefix=fastq_filename)
       # print("...done removing non-unique barcodes...")
    if maskOffByOne == True:
        masked_unique_split_loc_dict = MaskOffByOne(split_loc_dict, mapped_reads, outfileprefix=fastq_filename)
        #print("...done masking off-by-one barcodes...")

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
                for barcode in split_loc_dict[chrom][strand][pos]:
                    if split_loc_dict[chrom][strand][pos][barcode] > min_read_cutoff:
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

    wf = open(out_dir+fastq_filename+"_mapping_stats", 'a')
    wf.writelines("\n\nAnnotating insertions: \n")
    wf.writelines("reads mapping to NC plasmid backbone: "+str(plasmid_reads)+" ("+str(100*plasmid_reads/mapped_reads)+"%) of mapped reads\n")
    wf.writelines("reads mapping to TN right border: "+str(Rtn_reads)+" ("+str(100*Rtn_reads/mapped_reads)+"%) of mapped reads\n")
    wf.writelines("reads mapping to TN left border: "+str(Ltn_reads)+" ("+str(100*Ltn_reads/mapped_reads)+"%) of mapped reads\n")

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
    print("===merging_all_tnseq==")
    duplicate_reads = 0
    if len(loc_dicts)==1:
        return
    else:
        merged_dict = loc_dicts[0]
        duplicate_insert_dict = {}
        for loc_dict in loc_dicts[1:]:
            for chrom in loc_dict:
                if chrom not in merged_dict:
                    merged_dict[chrom] = {'+' : {}, '-' : {}} 
                for strand in loc_dict[chrom]:
                    for pos in loc_dict[chrom][strand]:
                        for barcode in loc_dict[chrom][strand][pos]:
                            total = loc_dict[chrom][strand][pos][barcode]
                            if pos not in merged_dict[chrom][strand]: # add an insert to the merged dict if not already one at that position
                                merged_dict[chrom][strand][pos]={}
                                merged_dict[chrom][strand][pos][barcode]=loc_dict[chrom][strand][pos][barcode]
                            elif barcode in merged_dict[chrom][strand][pos]: #if same barcode, add to total
                                merged_dict[chrom][strand][pos][barcode]+=total
                                duplicate_reads+=1
                                if barcode in duplicate_insert_dict:
                                    duplicate_insert_dict[barcode][1]+=1
                                else:
                                    duplicate_insert_dict[barcode]=[chrom+strand+str(pos), 2]
                            else:
                                merged_dict[chrom][strand][pos][barcode] = loc_dict[chrom][strand][pos][barcode]

                            
        
        wf=open(out_dir+merged_filename+"_duplicate_inserts", 'w')
        wf.writelines("combined "+str(duplicate_reads)+" duplicate inserts (same barcode at same position) from different fastq\n\n")
        wf.writelines("barcode\tposition\tnumber of times found\n")
        for dupins in duplicate_insert_dict:
            wf.writelines(dupins+"\t"+duplicate_insert_dict[dupins][0]+"\t"+str(duplicate_insert_dict[dupins][1])+"\n")
            
        wf.close()
        
        
        merged_dict = clean_barcodes(merged_dict, reads_remaining,single_barcode=single_barcode, maskOffByOne=maskOffByOne, outfileprefix=merged_filename)
        print("combined "+str(duplicate_reads)+" duplicate inserts (same barcode at same position) from different fastq")

        wf = open(out_dir+merged_filename+"_mapping_stats", 'a')
        wf.writelines('\n==='+merged_filename+'===\n')
        wf.writelines("BASH COMMAND\n"+command_line_string+ "\nRUN DATETIME\n" + now.strftime("%m/%d/%Y, %H:%M:%S")+"\n")
        wf.writelines("\nMerging libraries:\n")
        wf.writelines("combined "+str(duplicate_reads)+" duplicate inserts (same barcode at same position) from different fastq")
        wf.close()

        
        
        Annotate_insetions(merged_dict, reads_remaining, merged_filename) # identify insertions in genes, and write the final pooled outfile for a given fastq
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
    fastq_filename = read_file.split("/")[-1] #get file name

    print('\n==='+fastq_filename+'===\n')
    wf = open(out_dir+fastq_filename+"_mapping_stats", 'a')
    wf.writelines('\n==='+fastq_filename+'===\n')
    wf.writelines(command_line_string+ " run on " + now.strftime("%m/%d/%Y, %H:%M:%S")+"\n")
    wf.close()

    num_parsed_reads = Filter_Reads(read_file)  ## filters out reads that don't have tn sequence, writes the genomic portion of the remaining reads to a new file
    Map_Reads()  ## maps reads
    print("...done mapping...")

    loc_dict, reads_remaining = Filter_for_ID_and_multimapping_and_pool(num_parsed_reads)  # filters out reads below the identity threshold, that map to multiple locations and then pools insertions 

    loc_dict = Combine_near_mappings(loc_dict, reads_remaining) # combine insertions within 3bp of each other

    loc_dict = clean_barcodes(loc_dict, reads_remaining,single_barcode=single_barcode, maskOffByOne=maskOffByOne,outfileprefix=fastq_filename)
    Annotate_insetions(loc_dict, reads_remaining, fastq_filename) # identify insertions in genes, and write the final pooled outfile for a given fastq
    print ("...done annotating")
    loc_dicts.append(loc_dict)

if len(loc_dicts)>1:
    print()
    merged_filename = "merged_tnseq_"+now.strftime("%Y-%m-%d")+'.'
    merge_all_tnseq(loc_dicts,merged_filename) #merge dictionaries and make final pooled outfile for all fastq

end_time = time.time() #Time after it finished

print("Analysis took ", end_time-start_time, " seconds")

                
