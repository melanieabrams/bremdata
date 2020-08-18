import sys
from statistics import median
import random

'''
resample ranked h12 H2/H1 ratios by gene (accounting for gene size instead of essentiality)
 
'''
#Parameters

gff='/usr2/people/mabrams/Amended_Genomes/S288C/saccharomyces_cerevisiae_R64-1-1_20110208_withoutFasta.gff'
roman_numerals_in_gff=True

goi=['YLR397C','YGR098C', 'YMR168C','YKR054C','YDR180W',
     'YHR023W','YGR198W','YHR166C','YCR042C','YPL174C']

###
ranked_genes_G1 = sys.argv[1]

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
    gff_genes={}
    
    f = open(gfffile)
    lines=[]
    for line in f:
        if line[0]!='#' : #skip header rows
            row_data=line.split("\t")
            chrom=row_data[0]
            start=int(row_data[3])
            stop=int(row_data[4])
            gene_length=stop - start + 1
            info=row_data[8].split(";")
            yName=info[0].split('=')[1]
            if yName[0]=='Y' and len(yName)>5:
                gff_genes[yName]=gene_length
                
    f.close()
    
    return gff_genes


def test_goi_in_gff(goi):
    present=[]
    for gene in goi:
        if gene not in gff_genes:
            print(str(gene)+' not in gff')
        else:
            present.append(gene)
    print("tested the following goi: "+' '.join(present))
    return present

def binGFFbySize(gff_genes):
    binned={}
    for gene in gff_genes:
        rounded=round(gff_genes[gene],-3)
        if rounded in binned:
            binned[rounded].append(gene)
        else:
            binned[rounded]=[gene]

    #print(binned)
    return binned
    


def ParseRankedGenes(ranked_genes_g1):
    '''
    Input: ranked h12 file
    Output: dictionary of gene:G1
    '''

    ranked_dict={}      
    #add position and coverage for each base to a dictionary for that chromosome
    f = open(ranked_genes_g1)
    next(f)
    for line in f:
        row_data = line.strip().split("\t")
       # print(row_data)
        if row_data[0][0]=="Y" and len(row_data[0])>4:
            gene=row_data[0]
            G1=float(row_data[1])
            ranked_dict[gene]=G1
    f.close()

    for gene in goi:
        print("G1 of "+gene+":\t"+str(ranked_dict[gene]))

##    print(ranked_dict)
##    print(ranked_dict['YLR397C'])

    return ranked_dict


def resample(G1_dict, gff_genes, gff_bins, n_resample=10000):
    '''resamples from genes of same rounded length (to nearest kb)'''

    #get goi G1s and see if they're essential
    goi_G1s=[]
    for gene in goi:
        goi_G1s.append(G1_dict[gene])

    num_goi=len(goi_G1s)

    print('number of goi : '+str(num_goi))
        
    my_median=median(goi_G1s)  

    greater_medians=0.0
    random_medians=[]
        

    
    #add the correct number of genes from the corresponding size bin, and then see how many random groups >median
    for i in range(n_resample):
        random_G1s=[]
        for gene in goi:
            rounded=round(gff_genes[gene],-3)
            random_gene=random.choice(gff_bins[rounded])
            while random_gene not in G1_dict:
                random_gene=random.choice(gff_bins[rounded])
            random_G1s.append(G1_dict[random_gene])
        random_median=median(random_G1s)
        random_medians.append(random_median)
        if random_median>=my_median:
            greater_medians+=1

    rv=greater_medians/float(n_resample)
    print('rv: '+str(rv))
    print('n_resample: '+str(n_resample))
    print('median random median: ' +str(median(random_medians))+', goi median: '+str(my_median))
        
    
    



##RUN###

gff_genes=ParseFromGFF(gff)
tested_goi=test_goi_in_gff(goi)

gff_bins=binGFFbySize(gff_genes)

G1_dict=ParseRankedGenes(ranked_genes_G1)

resample(G1_dict,gff_genes,gff_bins)
