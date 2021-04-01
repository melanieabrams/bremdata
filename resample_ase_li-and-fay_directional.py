import sys
from statistics import median
import random


import matplotlib.pyplot as plt
import seaborn as sns


'''
resample ranked Supplementary Li and Fay 2017 outputs by gene,for directional divergence
(enrichment for more allele-and-species-specific upregulation across goi OR downregulation across goi than random samples)

 
'''
#Parameters
gff='saccharomyces_cerevisiae_R64-1-1_20110208_withoutFasta.gff'
roman_numerals_in_gff=True

essential_genes='essential.csv'


goi=['YLR397C','YGR098C','YMR168C','YKR054C',
        'YHR023W','YDR180W','YPL174C','YPR164W','YCR042C',
        'YMR016C','YJR135C','YJL025W','YDR443C',
        'YKL134C'] 

genes_ASE = 'Li_and_Fay_p37_allelefx.csv'

plot=True
savename='Li-and-Fay_directional_divergence'

###Functions
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

def ParseGenes(genes_ASE):
    '''
    Input: reformatted supplement of fay, with log2(avg sc reads 37/avg su reads 37) in column J
    Output: dictionary of gene:log2_sc_sp
    '''

    gene_dict={}      
    #add position and coverage for each base to a dictionary for that chromosome
    f = open(genes_ASE)
    next(f)
    for line in f:
        row_data = line.strip().split(",")
        gene=row_data[0]
        log2_sc_sp=float(row_data[1]) #cis-reg data
        gene_dict[gene]=log2_sc_sp
    f.close()

    for gene in goi:
        if gene not in gene_dict:
            gene_dict[gene]=None
        print("log2_sc_sp of "+gene+": "+str(gene_dict[gene]))

    return gene_dict


def resample(goi,log2_sc_sp_dict, essential_genes, n_resample=10000):
    '''resamples same num essential and nonessential genes'''

    #get goi log2_sc_sps and see if they're essential
    goi_log2_sc_sps=[]
    num_essential_goi=0
    for gene in goi:
        goi_log2_sc_sps.append(log2_sc_sp_dict[gene])
        if gene in essential_genes:
            num_essential_goi+=1
    num_goi=len(goi_log2_sc_sps)
    num_nonessential_goi=num_goi-num_essential_goi

    print('number essential goi : '+str(num_essential_goi))
    print('number nonessential goi: '+str(num_nonessential_goi))
        
    my_median=median(goi_log2_sc_sps)

    #split all tested log2_sc_sps into essential and nonessential
    ess_log2_sc_spdict={}
    noness_log2_sc_spdict={}
    for gene in log2_sc_sp_dict:
        if gene in essential_genes:
            ess_log2_sc_spdict[gene]=log2_sc_sp_dict[gene]
        else:
            noness_log2_sc_spdict[gene]=log2_sc_sp_dict[gene]

    

    greater_medians=0.0
    random_medians=[]

    #add the correct number of essential and nonessential random genes, and then see how many random groups >median
    for i in range(n_resample):
        random_log2_sc_sps=[]
        for n in range(num_essential_goi):
            random_gene=random.choice(list(ess_log2_sc_spdict))
            random_log2_sc_sps.append(log2_sc_sp_dict[random_gene])
        for n in range(num_nonessential_goi):
            random_gene=random.choice(list(noness_log2_sc_spdict))
            random_log2_sc_sps.append(log2_sc_sp_dict[random_gene])
        random_median=median(random_log2_sc_sps)
        random_medians.append(random_median)
        if abs(random_median)>=abs(my_median): #absolute values of medians - larger than OR smaller than random samples  
            greater_medians+=1

    rv=greater_medians/float(n_resample)
    print('rv: '+str(rv))
    print('n_resample: '+str(n_resample))
    print('median random median: ' +str(median(random_medians))+', goi median: '+str(my_median))

    if plot==True:  
        sns.displot(random_medians, kind='kde')
        plt.axvline(x=my_median, color='r', label = 'GOI median')
        plt.xlabel('Sample medians')
        plt.ylabel('Frequency')
        plt.title(savename+'\nresampling distribution, p= ' + str(rv))
        #plt.show()
        plt.tight_layout()
        plt.savefig(savename+'.png')
        
    
    



##RUN###

print("HS log2_sc/sp in cis\n")


gff_genes=ParseFromGFF(gff)
essential_genes=ParseEssential(essential_genes)
log2_sc_sp_dict=ParseGenes(genes_ASE)

tested_goi=test_goi_in_gff(goi)


resample(tested_goi,log2_sc_sp_dict,essential_genes)
