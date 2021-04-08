import sys
import subprocess as sp

# PARAMETERS # 


z1_dyn1 = "sgrp2_dyn1_z1.fasta"
d1373_dyn1 = "sgrp2_dyn1_d1373.fasta"
fastq_filename='IT016_S82_L007_R1_001.fastq'
blat_path='/usr2/people/mabrams/bin/blat/blat'


# BEGIN FUNCTIONS #

def Parse_Reads(fastq_file):  ## reformat so blat can use

    wf = open(fastq_filename+'_parsed_reads','w')  # outfile for the reads that will be mapped
  
    # counts for summary stats #

    tot_reads = 0.0  
    line_count = 0
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


                    
            # writes a new fastq file that has the read and header with the > # 

            wf.writelines(">"+header)
            wf.writelines(read+"\n")

            
    wf.close()    

  

    return 


def Map_Reads(mapping_genome,strain): # uses blat to map the reads
    
    cmd = [blat_path, mapping_genome, fastq_filename+"_parsed_reads", fastq_filename+"_mapped_reads_"+strain, "-minIdentity=95", "-tileSize=12", "-out=blast8" ]

    print(' '.join(cmd))
    
    sp.call(cmd)

                  

#### START PROGRAM ####
Parse_Reads(fastq_filename)  ## writes the reads from the fastq to a new file (getting rid of headers like ' @K00364:227:HGJ3LBBXY:7:1101:1458:1455')

Map_Reads(z1_dyn1,'Z1')  ## maps reads for spar
Map_Reads(d1373_dyn1,'D1373')  ## maps reads

print("done mapping")




