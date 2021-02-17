import pandas as pd
import numpy as np
import collections
import random

##Parameters###

go_file='go_terms.csv' #source: http://geneontology.org/docs/download-go-annotations/
essential_file='essential.csv' #source: http://www-sequence.stanford.edu/group/yeast_deletion_project/Essential_ORFs.txt

goi=['YGR198W','YLR397C','YBR136W','YGR140W','YNL049C','YKL017C','YGL082W',
     'YMR125W','YDR508C','YMR207C','YDR375C','YDR180W','YKL197C','YDR318W',
     'YMR094W','YOR326W','YBR081C','YPR049C','YIL152W','YER151C','YJR107W',
     'YAL026C','YDR456W','YLR141W','YPL268W','YDR235W']

#note: this go term enrichment script is modified from Dubin et al., 2020 (https://github.com/clairedubin/peroxisome_evolution/blob/68589c20c2e637d645ab4e79511db3ea992b121f/gene_expression/expression_analysis_artieri.ipynb)

###Functions###

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

def ParseEssential(essential_file):
    '''load list of essential genes as essential_genes'''
    print('...parsing essential genes...')
    essential = pd.read_csv(essential_file, header=None)
    essential[1] = essential[1].str.strip('\t')
    essential_genes = essential[1].str.strip(' ').tolist()
    return essential_genes



def resampleAllGo(go_term_groups,goi,go_terms,essential_count,non_essential_count,n=10000,save_intermediate=True):
    print('...resampling all GO terms...')

    
    #split go terms into essential and nonessential, and if in goi
    go_terms['essential']=go_terms.index.isin(essential_genes)
    go_terms['in_goi']=go_terms.index.isin(goi)
    
    essential_go_df=go_terms[go_terms['essential']==True]
    non_essential_go_df=go_terms[go_terms['essential']==False]   
    goi_go_df=go_terms[go_terms['in_goi']==True]

    print('- saving intermediate?: '+str(save_intermediate))
    if save_intermediate==True:
        goi_go_df.to_csv('goi_go_df.tsv', sep='\t', index=True)
    

    #get all go terms represented by more than one gene in goi
    all_goi_go_terms=goi_go_df.go_term.tolist()
    print('- number of go terms among goi (including dup.): '+str(len(all_goi_go_terms)))
    goi_go_terms=list(set(all_goi_go_terms))
    print('- number of go terms among goi (without dup.): '+str(len(goi_go_terms)))
    goi_go_dupes=[item for item, count in collections.Counter(all_goi_go_terms).items() if count > 1]
    print('-number of duplicated go terms among goi): '+str(len(goi_go_dupes)))
    #print('terms: \n',goi_go_dupes)
    ##count go terms in goi
    goi_go_counts=goi_go_df.go_term.value_counts()
    ##print(goi_go_counts)
 
    ##get n random samples of the same # of essential and nonessential genes
    print('- number of random samples: '+str(n))
    samples=[]
    for i in range(n):
        random_sample = list(random.sample(essential_go_df.index.values.tolist(),essential_count))
        random_sample+= list(random.sample(non_essential_go_df.index.values.tolist(),non_essential_count))
        go_terms['in_random']=go_terms.index.isin(random_sample)
        samples.append(go_terms[go_terms['in_random']==True])
    print(samples[0])
    go_terms=go_terms.drop(columns=['in_random'])
                           
    ##count go terms in each random sample
    sample_counts=[]
    for sample in samples:
        sample_counts.append(sample.go_term.value_counts())
    ##print('sample counts[:2]',sample_counts[:2])
        
    
    ##resample each go term duplicated in the goi with all n samples and the goi:
    rv_dict={}
    counter=0
    for term in goi_go_dupes:
        ##get number of that goi in go term

        n_go_goi=goi_go_counts[term]

        ##get number of genes in each random sample in go term, compare to goi
        n_samples_greater_or_equal_to_goi=0
        counts_from_random_samples=[]
        for sample_count in sample_counts:
            try:
                n_go_sample=sample_count[term]
            except KeyError:
                n_go_sample=0
##            if term=='GO:0005737':
##                print(n_go_goi,n_go_sample)
            if n_go_sample>=n_go_goi:
                n_samples_greater_or_equal_to_goi+=1
            counts_from_random_samples.append(n_go_sample)
        median_random=np.median(counts_from_random_samples)
            
        rv_dict[term]=[n_go_goi,median_random,float(n_samples_greater_or_equal_to_goi)/float(n)]
        


##    print(rv_dict['GO:0005737'])
    rv_df=pd.DataFrame.from_dict(rv_dict, orient='index',columns=['count_in_goi','median_count_in_random_sample','raw_rv'])
    rv_df['corrected_rv']=rv_df['raw_rv']*len(goi_go_dupes) #multiply rv by # of hypotheses
    rv_df.to_csv('resample_v2_goi_go_df.tsv', sep='\t', index=True)
##    print(rv_df)
  
    return
    

if __name__ == "__main__":

    #parse external data files
    go_terms=ParseGO(go_file)
    essential_genes=ParseEssential(essential_file)

    
    #split goi into essential and non essential genes

    essential_goi = [gene for gene in goi if gene in essential_genes]
    essential_count=len(essential_goi)
    print('\nessential goi ('+str(essential_count)+'):')
    print(essential_goi)
    non_essential_goi=[gene for gene in goi if gene not in essential_genes]
    non_essential_count=len(non_essential_goi)
    print('\nnonessential goi ('+str(non_essential_count)+'):')
    print(non_essential_goi)


    #add essentiality info to the the go term df
    print('\n...adding essentiality info for all genes, counting for go terms...')

    #print(go_terms)
    
    #count # essential genes in each go term
    go_term_groups = go_terms.groupby('go_term').count()
    go_term_groups.drop(columns=['sgd_name','gene_desc'],inplace=True)
    #print(go_term_groups)

    #resample
    resampleAllGo(go_term_groups,goi,go_terms,essential_count,non_essential_count,n=10000,save_intermediate=False)
    

    
