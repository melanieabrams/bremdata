import sys

# PARAMETERS # 
gff = sys.argv[1]
outputfilename = sys.argv[2]
WT_d_bedgraph =sys.argv[3]
WT_vcf = sys.argv[4]
mut_vcfs = sys.argv[5:]

if len(sys.argv)<5:
    print('USAGE: python3 gff outputfilename WT_d_bedgraph WT_vcf mut_vcf1 mutvcf2...')

# REF FILES #
sgdflatfile="/Users/Melanie/Desktop/Melanie/NGS/SGD_features.tab"
essentialorfs="/Users/Melanie/Desktop/Melanie/NGS/Essential_ORFs.txt"

# FUNCTIONS #

def ParseFromFile(vcf):
    '''
    Parses .vcf files to organize per chromosome, each mutation position and its metadata
    Input: .vcf 
    Output: {chrom#:[{pos1: [ref,alt,qual,mutType]}, {pos1: [ref,alt,qual,mutType]}...]}
    '''
    #build empty dictionary of chromosomes
    var_dict={}
    for i in range(0,16):
        var_dict[i+1]={}
         
    f = open(vcf)
    for line in f:
        if line[0]=="#":
            continue
        else:
            row_data =line.strip().split("\t")
            chrom = int(row_data[0][2:])
            pos = int(row_data[1])
##            if pos == 388357:
##                print("YDL039C SNP1 in "+vcf)
            ref = row_data[3]
            alt = row_data[4]
            qual = int(row_data[5])
            if row_data[7][0:5] == "INDEL":
                mutType = "INDEL"
            else:
                mutType = "SNP"
            var_dict[chrom][pos]=[ref,alt,pos,mutType]
    return var_dict

def ParseBedgraph(d_bedgraph):
    '''
    Input: -d bedgraph (i.e., with a value per
    Output: {'chrom':[[positions],[coverages]]
    '''
    #build empty dictionary of chromosomes
    coverage_dict={}
    for i in range(0,16):
        coverage_dict[i+1]=[],[]
         
    #add position and coverage for each base to a dictionary for that chromosome
    f = open(d_bedgraph)
    for line in f:
        row_data = line.strip().split("\t")
        chromosome = int(row_data[0][2:])
        pos = int(row_data[1])
        cov = int(row_data[2])
        coverage_dict[chromosome][0].append(pos)
        coverage_dict[chromosome][1].append(cov)
    #print(coverage_dict[1][1][0])
    return coverage_dict

def ParseSGDFlat(sgdflatfile):
    '''
    Parses SGD features flat file
    Input: SGD_features.tab file
    Output: 
    '''

    descr_dict={}
    
    f = open(sgdflatfile)
    for line in f:
        row_data =line.split("\t")
        yName=row_data[3]
        descr=row_data[15]
        descr_dict[yName]=descr
    return descr_dict

def ParseEssentials(essentialorfs):
    '''
    Parses essential genes yNames
    Input: SGD_features.tab file
    Output: 
    '''

    essentials=[]
    f = open(essentialorfs)
    for line in f:
        row_data =line.split("\t")
        if len(row_data)>1:
            yName= row_data[1].replace(" ","")
            essentials.append(yName)
    return essentials


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
    for i in range(0,16):
        hits[i+1] = []

    for chrom in vcf:
        for pos in vcf[chrom]:
 #           if pos == 388357:
 #               print 'YDL039 snp1'
            for gene in gff_dict[chrom]:
                start= gff_dict[chrom][gene][0]
                end = gff_dict[chrom][gene][1]
##                if start == 387158:
##                    if pos == 388357:
##                        print pos, start,end
##                        print pos>= start
##                        print pos <=end
##                        print gene
                if pos >= start and pos <= end:
##                    if pos == 388357:
##                        print pos, start,end
##                        print pos>= start
##                        print pos <=end
##                        print gene
##                        vcf_hit = vcf[chrom][pos]
##                        gff_hit = [1,gene]
##                        vcf_hit+=(gff_hit)
##                        hits[chrom].append(vcf_hit)
##                        print hits
                   # print pos, start, end,gene
                    vcf_hit = vcf[chrom][pos]
                    gff_hit = [gene]
                    vcf_hit+=(gff_hit)
                    hits[chrom].append(vcf_hit)

    return hits

def combine_vars(hitlist):
    all_hits = {}
    for i in range(0,16):
        all_hits[i+1] = {}

    for hit in hitlist:
        for chrom in hit:
            for var in hit[chrom]:
                ORF = var[-1]
                if ORF in all_hits[chrom]:
                    all_hits[chrom][ORF][0]+=1
                    all_hits[chrom][ORF].append(var[0:4])
                else:
                    all_hits[chrom][ORF] = [1,var[0:4]]

    return all_hits


def remove_WT_and_clean(all_hits, wt_var, threshold=3):
    clean_hits = {}
    for i in range(0,16):
        clean_hits[i+1] = {}
        
    for chrom in all_hits:
        for gene in all_hits[chrom]:
            if all_hits[chrom][gene][0] >= threshold:
                for var in wt_var[chrom]:
                    if gene in var:
                        continue
                    else:
                        clean_hits[chrom][gene]=all_hits[chrom][gene]
    return clean_hits

def remove_poor_coverage(clean_hits, wt_bg,covthreshold=10):
    clean_covered_hits = {}
    for i in range(0,16):
        clean_covered_hits[i+1] = {}
        
    covered_hits=clean_hits.copy()
    for chrom in covered_hits:
        for gene in covered_hits[chrom]:
            hitcount=covered_hits[chrom][gene][0]
            for var in covered_hits[chrom][gene][1:]:
                pos=var[2]
                if wt_bg[chrom][1][pos]<covthreshold:
                    covered_hits[chrom][gene].remove(var)
                    hitcount-=1

    for chrom in covered_hits:
        for gene in covered_hits[chrom]:
            if len(covered_hits[chrom][gene])>1:
                clean_covered_hits[chrom][gene]=covered_hits[chrom][gene]
            
    return clean_covered_hits
    

def print_combo(cleaned_combo):
    for chrom in cleaned_combo:
        print("Chromosome: "+str(chrom))
        for gene in cleaned_combo[chrom]:
            print("gene: "+gene)
            print(cleaned_combo[chrom][gene])
        print
    return None


def save_cleaned_var(cleaned_combo):
    SGD_dict=ParseSGDFlat(sgdflatfile)
    filename="cleaned_variants_"+outputfilename+".txt"
    f=open(filename,"w+")
    f.write ("gene" + "\t" + "hits" + "\t" +"type" + "\t" + "descr" +"\t" + "ref" +"\t"+"alt" +"\t"+"pos"+"\n")
    for chrom in cleaned_combo:
        for gene in cleaned_combo[chrom]:
            
            mutStr=""
            for variant in cleaned_combo[chrom][gene][1:]:
                mutStr+="\t"+variant[0]+"\t"+variant[1] + "\t" + str(chrom)+":"+ str(variant[2])
                               
            mutType= cleaned_combo[chrom][gene][1][-1] #assign as SNP, indel, or mix between the two
            for variant in cleaned_combo[chrom][gene][2:]:
                if variant[-1]!=mutType:
                    mutType = "MIX"

            if mutType == "SNP":
                mutType+="s"

            if mutType == "INDEL":  #assign as fs or non fs indel
                mods=[]
                for variant in cleaned_combo[chrom][gene][2:]:
                    mods.append(len(variant[0])%3.0 - len(variant[1])%3.0)
                for mod in mods:
                    if mod != 0.0:
                        if 0.0 in mods:
                            mutType = "INDELs(mixed)"
                        else:
                            mutType = "INDELs(fs)"
                    else:
                        mutType = "INDELs(no fs)"
                
            try:
                f.write(gene + "\t" + str(cleaned_combo[chrom][gene][0]) + "\t" +
                        mutType + "\t" + SGD_dict[gene].strip()+ mutStr +"\n")
            except KeyError:
                f.write(gene + "\t" + str(cleaned_combo[chrom][gene][0]) + "\t" +
                        mutType+mutStr+ "\n")
    f.close()

def save_essential_var(cleaned_combo):
    SGD_dict=ParseSGDFlat(sgdflatfile)
    essentialgenes=ParseEssentials(essentialorfs)
    filename="essential_cleaned_variants_"+outputfilename+".txt"
    f=open(filename,"w+")
    f.write ("gene" + "\t" + "hits" + "\t" +"type" + "\t" +"descr" + "\t"+"ref" +"\t"+"alt"+"\t"+"pos"+"\n")
    for chrom in cleaned_combo:
        for gene in cleaned_combo[chrom]:
            if gene in essentialgenes:
                mutStr=""
                for variant in cleaned_combo[chrom][gene][1:]:
                   mutStr+="\t"+variant[0]+"\t"+variant[1] + "\t" + str(chrom)+":"+ str(variant[2])
                   
                mutType= cleaned_combo[chrom][gene][1][-1]
                for variant in cleaned_combo[chrom][gene][2:]:
                    if variant[-1]!=mutType:
                        mutType = "MIX"

                if mutType == "SNP":
                    mutType+="s"

                if mutType == "INDEL":  #assign as fs or non fs indel
                    mods=[]
                    for variant in cleaned_combo[chrom][gene][2:]:
                        mods.append(len(variant[0])%3.0 - len(variant[1])%3.0)
                    for mod in mods:
                        if mod != 0.0:
                            if 0.0 in mods:
                                mutType = "INDELs(mixed)"
                            else:
                                mutType = "INDELs(fs)"
                        else:
                            mutType = "INDELs(no fs)"
                        
                try:
                    f.write(gene + "\t" + str(cleaned_combo[chrom][gene][0]) + "\t" +
                            mutType + "\t" + SGD_dict[gene].strip()+ mutStr +"\n")
                except KeyError:
                    f.write(gene + "\t" + str(cleaned_combo[chrom][gene][0]) + "\t" +
                            mutType + mutStr +"\n")
    f.close()
#### START PROGRAM ####

#import vcfs
wild_var=ParseFromFile(WT_vcf)

mut_vars=[]
for mut_vcf in mut_vcfs:
    mut_vars.append(ParseFromFile(mut_vcf))

#import d_bedgraph
wt_bg = ParseBedgraph(WT_d_bedgraph)

print("done importing")

#build annotation dictionary
gff_dict=build_gff_dict(gff)

#create annotated vcfs
annotated_wt=Annotate_CDS_variants(gff_dict,wild_var)
#print annotated_wt

annotated_muts=[]
for mut_var in mut_vars:
    annotated_muts.append(Annotate_CDS_variants(gff_dict,mut_var))

#print annotated_muts[0]

print("done annotating")

#combine 
combined_mut_an=combine_vars(annotated_muts)

cleaned_combo=remove_WT_and_clean(combined_mut_an, annotated_wt, threshold=1)

cleaned_combo=remove_poor_coverage(cleaned_combo, wt_bg, covthreshold=20)


print("done comparing")
#print_combo(cleaned_combo) #raw print of output, grouped by chromosome, [hits, [variant1, qual, type],[variant2...]
save_cleaned_var(cleaned_combo) #saved txt file of all output
save_essential_var(cleaned_combo) #saved txt file of output corresponding to essential genes


####
##
##


