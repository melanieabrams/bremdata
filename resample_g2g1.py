import sys
from statistics import median
import random

'''
resample ranked G2/G1 (or H2/H1) ratios by gene (accounting for essentiality)
 
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

def ParseRankedGenes(ranked_genes_h12,goi,printGeneG2G1=True):
    '''
    Input: ranked h12 file
    Output: dictionary of gene:H2H1
    '''

    ranked_dict={}      
    #add position and coverage for each base to a dictionary for that chromosome
    f = open(ranked_genes_h12)
    next(f)
    for line in f:
        row_data = line.strip().split("\t")
       # print(row_data)
        if row_data[0][0]=="Y" and len(row_data[0])>4:
            gene=row_data[0]
            H12=float(row_data[1])
            H2H1=float(row_data[2])
            ranked_dict[gene]=H2H1
    f.close()

    for gene in goi:
        if printGeneG2G1==True:
            print("G2G1 of "+gene+": \t"+str(ranked_dict[gene]))

##    print(ranked_dict)
##    print(ranked_dict['YLR397C'])

    return ranked_dict


def resample(H2H1_dict, essential_genes, goi, n_resample=10000, printResampling=True):
    '''resamples same num essential and nonessential genes'''

    #get goi G2G1s and see if they're essential
    goi_H2H1s=[]
    num_essential_goi=0
    for gene in goi:
        goi_H2H1s.append(H2H1_dict[gene])
        if gene in essential_genes:
            num_essential_goi+=1
    num_goi=len(goi_H2H1s)
    num_nonessential_goi=num_goi-num_essential_goi

    if printResampling==True:
        print('number essential goi : \t'+str(num_essential_goi))
        print('number nonessential goi: \t'+str(num_nonessential_goi))
        
    my_median=median(goi_H2H1s)

    #split all tested H2H1s into essential and nonessential
    ess_H2H1dict={}
    noness_H2H1dict={}
    for gene in H2H1_dict:
        if gene in essential_genes:
            ess_H2H1dict[gene]=H2H1_dict[gene]
        else:
            noness_H2H1dict[gene]=H2H1_dict[gene]

    

    greater_medians=0.0
    random_medians=[]

    #add the correct number of essential and nonessential random genes, and then see how many random groups >median
    for i in range(n_resample):
        random_H2H1s=[]
        for n in range(num_essential_goi):
            random_gene=random.choice(list(ess_H2H1dict))
            random_H2H1s.append(H2H1_dict[random_gene])
        for n in range(num_nonessential_goi):
            random_gene=random.choice(list(noness_H2H1dict))
            random_H2H1s.append(H2H1_dict[random_gene])
        random_median=median(random_H2H1s)
        random_medians.append(random_median)
        if random_median>my_median:
            greater_medians+=1

    rv=greater_medians/float(n_resample)
    if printResampling==True:
        print('rv: '+str(rv))
        print('n_resample: '+str(n_resample))
        print('median random median: \t' +str(median(random_medians))+', goi median: \t'+str(my_median))
            
    return rv
    



##RUN###

if __name__ == "__main__":

    ranked_genes_h12 = sys.argv[1]

    gff_genes=ParseFromGFF(gff)
    tested_goi=test_goi_in_gff(goi)

    essential_genes=ParseEssential(essential_genes)
    H2H1_dict=ParseRankedGenes(ranked_genes_h12)
    resample(H2H1_dict,essential_genes,goi)

