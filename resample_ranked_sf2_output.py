import sys
from statistics import median
import random

'''
resample sweepfinder 2 outputs by gene
 
'''
#Parameters
gff='/usr2/people/mabrams/Amended_Genomes/Z1/Z1.gff'
#gff='Z1.gff'
save_prefix='Bergstrom2014_ranked_genes_resample'
goi=['YLR397C','YGR098C', 'YMR168C','YKR054C','YDR180W',
     'YHR023W','YGR198W','YHR166C','YCR042C','YPL174C']

###
ranked_genes_sf2 = sys.argv[1]


def ParseFromGFF(gfffile):
    '''
    Parses SGD features flat file
    Input: SGD_features.tab file
    Output: dict of {gene:[chrom,start,stop]}}
    '''

    gff_genes=[]
    
    f = open(gfffile)
    lines=[]
    for line in f:
        row_data=line.split("\t")
        chrom=row_data[0]
        start=int(row_data[3])
        stop=int(row_data[4])
        info=row_data[8].split(";")
        yName=info[0].split('=')[1]
        if yName[0]=='Y':
            gff_genes.append(yName)
            #ann_dict[yName]=[chrom,start,stop]
                
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

def ParseRankedGenes(ranked_genes_sf2):
    '''
    Input: ranked sf2 file
    Output: 
    '''

    ranked_dict={}      
    #add position and coverage for each base to a dictionary for that chromosome
    f = open(ranked_genes_sf2)
    next(f)
    for line in f:
        row_data = line.strip().split("\t")
       # print(row_data)
        if row_data[0][0]=="Y" and len(row_data[0])>4:
            gene=row_data[0]
            LR=float(row_data[1])
            ranked_dict[gene]=LR
    f.close()

##    print(ranked_dict)
##    print(ranked_dict['YLR397C'])

    return ranked_dict


def resample(LR_dict, n_resample=10000):
    goi_LRs=[]
    for gene in goi:
        goi_LRs.append(LR_dict[gene])
    my_median=median(goi_LRs)

    num_goi=len(goi_LRs)

    greater_medians=0.0

    for i in range(n_resample):
        random_LRs=[]
        for n in range(num_goi):
            random_gene=random.choice(list(LR_dict))
            random_LRs.append(LR_dict[random_gene])
        random_median=median(random_LRs)
        if random_median>my_median:
            greater_medians+=1

    rv=greater_medians/float(n_resample)
    print('rv: '+str(rv))
    print('n_resample: '+str(n_resample))
        
    
    



##RUN###

gff_genes=ParseFromGFF(gff)
tested_goi=test_goi_in_gff(goi)
LR_dict=ParseRankedGenes(ranked_genes_sf2)
resample(LR_dict)
