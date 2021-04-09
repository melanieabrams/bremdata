import pandas as pd
import numpy as np
import collections
import random
import statsmodels.stats.multitest as smm

import matplotlib.pyplot as plt
import seaborn as sns

pd.options.mode.chained_assignment = None  # default='warn'

##Parameters###  ##NEED TO FINISH THIS - GET LIST OF ALL GENES TESTED FOR ESSENTIALITY RESMAPLING

essential_file='essential.csv' #source: http://www-sequence.stanford.edu/group/yeast_deletion_project/Essential_ORFs.txt
gff='saccharomyces_cerevisiae_R64-1-1_20110208_withoutFasta.gff'
go_file='go_terms.csv' #source: http://geneontology.org/docs/download-go-annotations/


goi=['YGR198W','YLR397C','YBR136W','YGR140W','YNL049C','YKL017C','YGL082W',
     'YMR125W','YDR508C','YMR207C','YDR375C','YDR180W','YKL197C','YDR318W',
     'YMR094W','YOR326W','YBR081C','YPR049C','YIL152W','YER151C','YJR107W',
     'YAL026C','YDR456W','YLR141W','YPL268W','YDR235W']

plot=True #for histogram of rv
savename='resample_essential_excluding_cellcycle'

processes_to_exclude='GO:0007049'

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

def ParseGO(go_file):
    '''load go term data as go_terms, remove 3 broad go terms'''
    print('...parsing go terms...')
    go_terms = pd.read_csv(go_file, header=None)
    go_terms = go_terms.drop(columns=[0, 1,3,5,6,7,8,11,12,13,14,15])
    go_terms = go_terms.rename(columns={2: 'sgd_name', 4:'go_term', 9:'gene_desc', 10:'gene'})
    go_terms['gene'] = [i[0] for i in go_terms['gene'].str.split('|')]
    go_terms = go_terms.drop_duplicates()
    go_terms = go_terms[~go_terms['go_term'].isin(['GO:0005575', 'GO:0008150', 'GO:0003674'])]
    go_terms = go_terms.set_index('gene')
    #print(go_terms.head())
    return go_terms


def excludeProcess(processes_to_exclude,go_terms):
    '''exclude genes with a particular go term from analysis'''
    #get a list of genes with the excluded process
    exclude_go=go_terms[go_terms['go_term']==processes_to_exclude]
    exclude_gene_list=exclude_go.index.to_list()
    #exclude those from the gprocesses_to_excludeo_terms df                  
    filtered_go_terms=go_terms[~go_terms.index.isin(exclude_gene_list)]
    return filtered_go_terms,exclude_gene_list



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

    if plot==True:  
        sns.displot(random_samples_num_essential, kind='kde')
        plt.axvline(x=num_essential_goi, color='r', label = 'GOI median')
        plt.xlabel('Sample counts')
        plt.ylabel('Frequency')
        plt.title(savename+'\nresampling distribution, p= ' + str(rv/n))
        #plt.show()
        plt.tight_layout()
        plt.savefig(savename+'.png')
        
            
  
    return
    

if __name__ == "__main__":


    
    np.random.seed(777)
    random.seed(777)

    #parse external data files
    gff_genes=ParseFromGFF(gff)
    essential_genes=ParseEssential(essential_file)
    go_terms=ParseGO(go_file)

    #exclude cell cycle go process
    
    go_terms,exclude_gene_list=excludeProcess(processes_to_exclude,go_terms)
    goi=[gene for gene in goi if gene not in exclude_gene_list]
    print('goi with cell cycle genes excluded:')
    print(goi)

    #resample
    resampleEssential(goi,gff_genes,essential_genes,n=10000)

    
