import sys

# PARAMETERS # 

gff = sys.argv[1]
WT_vcf = sys.argv[2]
mut_vcfs = sys.argv[3:]

# FUNCTIONS #

def ParseFromFile(vcf):
    '''
    Parses .vcf files to organize per chromosome, each mutation position and its metadata
    Input: .vcf 
    Output: {chrom#:[{pos1: [ref,alt,qual,mutType]}, {pos1: [ref,alt,qual,mutType]}...]}
    '''
    #build empty dictionary of chromosomes
    var_dict={}
    for i in xrange(0,16):
        var_dict[i+1]=[]
         
    f = open(vcf)
    for line in f:
        if line[0]=="#":
            continue
        else:
            row_data =line.strip().split("\t")
            chrom = int(row_data[0][2:])
            pos = int(row_data[1])
            ref = row_data[3]
            alt = row_data[4]
            qual = int(row_data[5])
            if row_data[7][0:5] == "INDEL":
                mutType = "INDEL"
            else:
                mutType = "SNP"
            var_dict[chrom].append({pos:[ref,alt,qual,mutType]})
    return var_dict


def build_gff_dict(gff):
    '''
    Input: gff annotation file
    Output: {chrom#:{[ORF]:[start,end,strand]}}
    '''
    gff_dict = {}
    f = open(gff)

    for line in f:
        line=line.strip().split("\t")
        if line[2] == "CDS":
            chrom = int(line[0])
            start = int(line[3])
            end = int(line[4])
            strand = line[6]
            ORF = line[8].split(";")[4][4:]
            if ORF == '':
                ORF = "No SGD_"+str(start)+"_to_"+str(end)

        # put the annotation information in a dictionary # 
            if chrom not in gff_dict:
                gff_dict[chrom] = {}

            gff_dict[chrom][ORF] = [start,end,strand]

    return gff_dict

    
def Annotate_CDS_variants(gff_dict,vcf):
    hits = {}
    for i in xrange(0,16):
        hits[i+1] = []

    for chrom in vcf:
        for pos in vcf[chrom][0]:
            for gene in gff_dict[chrom]:
                start= gff_dict[chrom][gene][0]
                end = gff_dict[chrom][gene][1]
                if pos >= start and pos <= end:
                    vcf_hit = vcf[chrom][0][pos]
                    gff_hit = [1,gene]
                    vcf_hit+=(gff_hit)
                    hits[chrom].append(vcf_hit)
   
    return hits

def combine_vars(hitlist):
    all_hits = {}
    for i in xrange(0,16):
        all_hits[i+1] = {}

    for hit in hitlist:
        for chrom in hit:
            if len(hit[chrom])>0:
                ORF = hit[chrom][0][5]
                if ORF in all_hits[chrom]:
                    all_hits[chrom][ORF][0]+=1
                    all_hits[chrom][ORF].append(hit[chrom][0][0:4])
                else:
                    all_hits[chrom][ORF] = [1,hit[chrom][0][0:4]]

    return all_hits


def remove_WT_vars(all_hits, wt_var):
    clean_hits = {}
    for i in xrange(0,16):
        clean_hits[i+1] = {}
    for chrom in all_hits:
        for gene in all_hits[chrom]:
            for var in wt_var[chrom]:
                if gene in var:
                    continue
                else:
                    clean_hits[chrom][gene]=all_hits[chrom][gene]
    return all_hits

#### START PROGRAM ####

#import vcfs
wild_var=ParseFromFile(WT_vcf)

mut_vars=[]
for mut_vcf in mut_vcfs:
    mut_vars.append(ParseFromFile(mut_vcf))

print "done importing"

#build annotation dictionary
gff_dict=build_gff_dict(gff)

#create annotated vcfs
annotated_wt=Annotate_CDS_variants(gff_dict,wild_var)

annotated_muts=[]
for mut_var in mut_vars:
    annotated_muts.append(Annotate_CDS_variants(gff_dict,mut_var))

print "done annotating"

#combine 
combined_mut_an=combine_vars(annotated_muts)

combo_minus_WT=remove_WT_vars(combined_mut_an, annotated_wt)

print "done comparing"

for chrom in combo_minus_WT:
    print "Chromosome: "+str(chrom)
    for gene in combo_minus_WT[chrom]:
        print "gene: "+gene
        print combo_minus_WT[chrom][gene]
    print






