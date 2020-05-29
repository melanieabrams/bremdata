import sys
from statistics import median
import random

'''
resample ranked sweepfinder 2 outputs by gene
 
'''
#Parameters
#gff='/usr2/people/mabrams/Amended_Genomes/Z1/Z1.gff'
#gff='Z1.gff'
#gff='/usr2/people/mabrams/Amended_Genomes/D1373/DBVPG1373.gff'
#gff='DBVPG1373.gff'
gff='/usr2/people/mabrams/Amended_Genomes/S288C/saccharomyces_cerevisiae_R64-1-1_20110208_withoutFasta.gff'

essential_genes='/usr2/people/mabrams/Amended_Genomes/essential.csv'

goi=['YLR397C','YGR098C', 'YMR168C','YKR054C','YDR180W',
     'YHR023W','YCR042C','YNL172W']

###
ranked_genes_sf2 = sys.argv[1]


def ParseFromGFF(gfffile):
    '''
    Parses SGD features flat file
    Input: SGD_features.tab file
    Output: dict of {chrom:{gene:[start,stop]}}
    '''
    gff_genes=[]
    
    f = open(gfffile)
    lines=[]
    for line in f:
        if line[0]!='#': #skip header rows
            row_data=line.split('\t')
            chrom=row_data[0]
            start=int(row_data[3])
            stop=int(row_data[4])
            info=row_data[8].split(";")
            yName=info[0].split('=')[1]
            #print(yName)
            if yName[0]=='Y' and len(yName)>5:
                gff_genes.append(yName)
    f.close()
    
    return gff_genes

def ParseEssential(essentialfile):
    '''
    Parses csv where Ynames of essential genes is second column
    Output: list of genes
    '''

    ess_genes=[]
    
    f = open(essentialfile)
    lines=[]
    bad_chars=['\t',' ','"']
    for line in f:
        row_data=line.split(",")
        yName=row_data[1]
        for char in bad_chars:
            yName=yName.replace(char, '')
        ess_genes.append(yName)
                
    f.close()

    #print(ess_genes)
    return ess_genes

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
    Output: dictionary of gene:LR
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

    for gene in goi:
        print("LR of "+gene+": "+str(ranked_dict[gene]))

##    print(ranked_dict)
##    print(ranked_dict['YLR397C'])

    return ranked_dict


def resample(LR_dict, essential_genes, n_resample=10000):
    '''resamples same num essential and nonessential genes'''

    #get goi LRs and see if they're essential
    goi_LRs=[]
    num_essential_goi=0
    for gene in goi:
        goi_LRs.append(LR_dict[gene])
        if gene in essential_genes:
            num_essential_goi+=1
    num_goi=len(goi_LRs)
    num_nonessential_goi=num_goi-num_essential_goi

    print('number essential goi : '+str(num_essential_goi))
    print('number nonessential goi: '+str(num_nonessential_goi))
        
    my_median=median(goi_LRs)

    #split all tested LRs into essential and nonessential
    ess_LRdict={}
    noness_LRdict={}
    for gene in LR_dict:
        if gene in essential_genes:
            ess_LRdict[gene]=LR_dict[gene]
        else:
            noness_LRdict[gene]=LR_dict[gene]

    

    greater_medians=0.0
    random_medians=[]

    #add the correct number of essential and nonessential random genes, and then see how many random groups >median
    for i in range(n_resample):
        random_LRs=[]
        for n in range(num_essential_goi):
            random_gene=random.choice(list(ess_LRdict))
            random_LRs.append(LR_dict[random_gene])
        for n in range(num_nonessential_goi):
            random_gene=random.choice(list(noness_LRdict))
            random_LRs.append(LR_dict[random_gene])
        random_median=median(random_LRs)
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

essential_genes=ParseEssential(essential_genes)
LR_dict=ParseRankedGenes(ranked_genes_sf2)
resample(LR_dict,essential_genes)
