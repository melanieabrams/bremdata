import pandas as pd
import numpy as np
import collections
import random
import statsmodels.stats.multitest as smm

pd.options.mode.chained_assignment = None  # default='warn'

##Parameters###  ##NEED TO FINISH THIS - GET LIST OF ALL GENES TESTED FOR ESSENTIALITY RESMAPLING

essential_file='essential.csv' #source: http://www-sequence.stanford.edu/group/yeast_deletion_project/Essential_ORFs.txt
gff='saccharomyces_cerevisiae_R64-1-1_20110208_withoutFasta.gff'


goi=['YGR198W','YLR397C','YBR136W','YGR140W','YNL049C','YKL017C','YGL082W',
     'YMR125W','YDR508C','YMR207C','YDR375C','YDR180W','YKL197C','YDR318W',
     'YMR094W','YOR326W','YBR081C','YPR049C','YIL152W','YER151C','YJR107W',
     'YAL026C','YDR456W','YLR141W','YPL268W','YDR235W']

carlygoi=['YLR397C','YGR098C','YKR054C',
     'YHR023W','YDR180W','YPL174C','YCR042C',
     'YMR016C','YJR135C','YJL025W','YDR443C',
     'YKL134C'] #from popgen MS, excluding YMR168C b/c only one br

goi=goi+carlygoi


###Functions###

def ParseFromGFF(gfffile):
    '''
    Parses gff
    Output: list of genes
    '''

    gene_list=[]
    
    f = open(gfffile)
    lines=[]
    for line in f:
        if line[0]!='#': #skip header rows
            row_data=line.split("\t")
            info=row_data[8].split(";")
            yName=info[0].split('=')[1]
            #get complient yNames
            if yName[0]=='Y' and len(yName)>5:
                if yName[-1]=="W" or yName[-1]=="C":
                    gene_list.append(yName)
    f.close()
    
    return gene_list

def ParseEssential(essential_file):
    '''load list of essential genes as essential_genes'''
    print('...parsing essential genes...')
    essential = pd.read_csv(essential_file, header=None)
    essential[1] = essential[1].str.strip('\t')
    essential_genes = essential[1].str.strip(' ').tolist()
    return essential_genes


def CountEssential(gene_list,essential_genes):
    num_essential=len([gene for gene in gene_list if gene in essential_genes])
    return num_essential

def resampleEssential(goi,gff_genes,essential_genes,n=10000):
    #get goi essential
    num_goi=len(goi)
    num_essential_goi=CountEssential(goi,essential_genes)
    print(str(num_essential_goi)+' of '+str(num_goi)+' goi are essential')

    ##get n random samples of the same # genes as in the goi, from the gff genes, count # essential in them
    print('- number of random samples: '+str(n))
    random_samples_num_essential=[]
    for i in range(n):
        random_sample = list(random.sample(gff_genes,num_goi))
        random_essential=CountEssential(random_sample,essential_genes)
        random_samples_num_essential.append(random_essential)


    rv=0.0
    for random_essential in random_samples_num_essential:
        if random_essential>=num_essential_goi:
            rv+=1
    
    print('resampling value: '+str(rv/n))
            
  
    return
    

if __name__ == "__main__":

    
    np.random.seed(777)

    #parse external data files
    gff_genes=ParseFromGFF(gff)
    essential_genes=ParseEssential(essential_file)

    #resample
    resampleEssential(goi,gff_genes,essential_genes,n=10000)

    
