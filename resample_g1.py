import sys
from statistics import median
import random

'''
resample ranked g1 by gene (or H1) (accounting for essentiality)
 
'''
#Parameters

gff='/usr2/people/mabrams/Amended_Genomes/S288C/saccharomyces_cerevisiae_R64-1-1_20110208_withoutFasta.gff'

essential_genes='/usr2/people/mabrams/Amended_Genomes/essential.csv'


goi=['YLR397C','YGR098C','YMR168C','YKR054C',
    'YHR023W','YDR180W','YPL174C','YCR042C',
    'YMR016C','YJR135C','YJL025W','YDR443C',
    'YKL134C']
###


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
        if line[0]!='#' : #skip header rows
            row_data=line.split("\t")
            chrom=row_data[0]
            start=int(row_data[3])
            stop=int(row_data[4])
            info=row_data[8].split(";")
            yName=info[0].split('=')[1]
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

def test_goi_in_gff(goi,gff_genes,printGOI=True):
    present=[]
    for gene in goi:
        if gene not in gff_genes:
            print(str(gene)+' not in gff')
        else:
            present.append(gene)
    if printGOI==True:
        print("tested the following goi: "+' '.join(present))
    return present

def ParseRankedGenes(ranked_genes_g1,goi,printGeneG1=True):
    '''
    Input: ranked G1 file from .h12_h2h1
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
        if printGeneG1==True:
            print("G1 of "+gene+": \t"+str(ranked_dict[gene]))

##    print(ranked_dict)
##    print(ranked_dict['YLR397C'])

    return ranked_dict


def resample(G1_dict, essential_genes, goi, n_resample=10000,printResampling=True):
    '''resamples same num essential and nonessential genes'''

    #get goi G1s and see if they're essential
    goi_G1s=[]
    num_essential_goi=0
    for gene in goi:
        goi_G1s.append(G1_dict[gene])
        if gene in essential_genes:
            num_essential_goi+=1
    num_goi=len(goi_G1s)
    num_nonessential_goi=num_goi-num_essential_goi

    if printResampling==True:
        print('number essential goi : \t'+str(num_essential_goi))
        print('number nonessential goi: \t'+str(num_nonessential_goi))
        
    my_median=median(goi_G1s)

    #split all tested G1s into essential and nonessential
    ess_G1dict={}
    noness_G1dict={}
    for gene in G1_dict:
        if gene in essential_genes:
            ess_G1dict[gene]=G1_dict[gene]
        else:
            noness_G1dict[gene]=G1_dict[gene]

    

    greater_medians=0.0
    random_medians=[]

    #add the correct number of essential and nonessential random genes, and then see how many random groups >median
    for i in range(n_resample):
        random_G1s=[]
        for n in range(num_essential_goi):
            random_gene=random.choice(list(ess_G1dict))
            random_G1s.append(G1_dict[random_gene])
        for n in range(num_nonessential_goi):
            random_gene=random.choice(list(noness_G1dict))
            random_G1s.append(G1_dict[random_gene])
        random_median=median(random_G1s)
        random_medians.append(random_median)
        if random_median>=my_median:
            greater_medians+=1

    rv=greater_medians/float(n_resample)
    if printResampling==True:
        print('rv: '+str(rv))
        print('n_resample: '+str(n_resample))
        print('median random median: \t' +str(median(random_medians))+'\n, goi median: \t'+str(my_median))
            
    return rv
    



##RUN###

if __name__ == "__main__":

    
    ranked_genes_G1 = sys.argv[1]

    gff_genes=ParseFromGFF(gff)
    tested_goi=test_goi_in_gff(goi,gff_genes)

    essential_genes=ParseEssential(essential_genes)
    G1_dict=ParseRankedGenes(ranked_genes_G1,goi)
    rv=resample(G1_dict,essential_genes,goi)
