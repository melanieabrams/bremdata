import pandas as pd
import numpy as np
import collections
import random
import statsmodels.stats.multitest as smm

import networkx
import obonet

##Parameters###

go_file='go_terms.csv' #source: http://geneontology.org/docs/download-go-annotations/
essential_file='essential.csv' #source: http://www-sequence.stanford.edu/group/yeast_deletion_project/Essential_ORFs.txt
fxsize_file='37C_2.0_1.1_coeffvarANDreadcutoff_uncollapsedMBA_unpaired_ratiocvfilt_10.0_3.0.all_fxsizes'
sep='\t'

n=10000

min_go_count=5 #discard GO terms with fewer genes than this
max_go_count=200 #discard GO terms with more genes than this

low_n_file=None #default - have not parsed from low n file
low_n_file='resample_all_go_by_absolute_fxsize_n100.tsv' #if have already run the code with low n (n=100, then 1000), to speed up computation of high n p-values for candidates with rv<0.1

#running on thar as of 2:25PM on 5/4/21
#note: this go term enrichment script is modified from Dubin et al., 2020 (https://github.com/clairedubin/peroxisome_evolution/blob/68589c20c2e637d645ab4e79511db3ea992b121f/gene_expression/expression_analysis_artieri.ipynb)

###Functions###

def ParseFxSize(fxsize_file):
    '''return absolute effect sizes from barseq data'''
    fx_df=pd.read_csv(fxsize_file,sep=sep)
    fx_df['absolute_fx']=abs(fx_df['sc_mean']-fx_df['sp_mean'])
    keep_col=['gene','absolute_fx']
    fx_df=fx_df.drop(columns=[col for col in fx_df.columns if col not in keep_col])
    return fx_df
    

def ParseGO(go_file):
    '''load go term data as go_terms, remove 3 broad go terms'''
    print('...parsing go terms...')
    go_terms = pd.read_csv(go_file, header=None)
    go_terms = go_terms.drop(columns=[0, 1,3,5,6,7,8,11,12,13,14,15])
    go_terms = go_terms.rename(columns={2: 'sgd_name', 4:'go_term', 9:'gene_desc', 10:'gene'})
    go_terms['gene'] = [i[0] for i in go_terms['gene'].str.split('|')]
    go_terms = go_terms.drop_duplicates()
    go_terms = go_terms[~go_terms['go_term'].isin(['GO:0005575', 'GO:0008150', 'GO:0003674'])] 
    #go_terms = go_terms.set_index('gene')
    #print(go_terms.head())
    return go_terms

def ParseEssential(essential_file):
    '''load list of essential genes as essential_genes'''
    print('...parsing essential genes...')
    essential = pd.read_csv(essential_file, header=None)
    essential[1] = essential[1].str.strip('\t')
    essential_genes = essential[1].str.strip(' ').tolist()
    return essential_genes

def resampleOneGo(go,go_terms,essential_go_counts,non_essential_go_counts,fx_df,n=10000):
    #count number of essential and nonessential genes in this go term
    #print(essential_go_counts)
    try:
        this_go_essential_count=essential_go_counts._get_value(go,'count')
    except KeyError:
        this_go_essential_count=0
    try:
        this_go_nonessential_count=non_essential_go_counts._get_value(go,'count')
    except KeyError:
        this_go_nonessential_count=0

    #this_go_nonessential_count=len(set(non_essential_go_df[non_essential_go_df['go_term'].isin(genes_in_this_go)]))

   
    #get median fitness of this go term
    genes_in_this_go=go_terms[go_terms['go_term']==go]['gene'].to_list()
    this_fx_df=fx_df[fx_df['gene'].isin(genes_in_this_go)]
    #print(this_fx_df)
    this_median_fx=this_fx_df['absolute_fx'].median()


    #get random samples, tabulate their median fitness, and add to rv if greater
    rv=0.0
    for i in range(n):
        random_sample = list(random.sample(essential_go_df['gene'].tolist(),this_go_essential_count))
        random_sample+= list(random.sample(non_essential_go_df['gene'].tolist(),this_go_nonessential_count))
        #print(random_sample)
        random_fx_df=fx_df[fx_df['gene'].isin(random_sample)]
        #print(this_go_essential_count,this_go_nonessential_count)
        #print(random_fx_df)
        random_median_fx=random_fx_df['absolute_fx'].median()
        #print(this_median_fx,random_median_fx)
        if random_median_fx > this_median_fx:
            rv+=1
    
    return rv,this_median_fx,random_median_fx,this_go_essential_count,this_go_nonessential_count



def resampleAllGo(go_terms,fx_df,essential_go_counts,non_essential_go_counts,low_n_file=low_n_file,n=10000):
    print('...resampling all GO terms...')

    each_go=list(set(go_terms['go_term']))

    if low_n_file!=None: #if already have rv calculated from a file with lower resampliing
        low_n_df=pd.read_csv(low_n_file,sep=sep)
        low_n_df.drop(columns=['bh_rv'],inplace=True)#will be recalculating multiple hypothesis with higher n rv
        go_for_further_resampling=low_n_df[low_n_df['raw_rv']<=0.1]['Unnamed: 0'].to_list()
        print('there are ' +str(len(go_for_further_resampling))+' go terms with rv <0.1 for closer resampling')
        low_n_df.set_index('Unnamed: 0',inplace=True)
        #low_n_df=low_n_df.head()
        

        rv_dict=low_n_df.T.to_dict('list')

        counter=0
        for go in go_for_further_resampling:
            counter+=1
            if counter%100==0:
                print(counter)
            rv,this_median_fx,random_median_fx,num_essential,num_nonessential=resampleOneGo(go,go_terms,essential_go_counts,non_essential_go_counts,fx_df,n=n)
            rv_dict[go]=rv/n,this_median_fx,random_median_fx,num_essential,num_nonessential
    

    else:

        rv_dict={}
        #print(each_go)

        counter=0
        for go in each_go:
            
            counter+=1
            if counter%100==0:
                print(counter)
            rv,this_median_fx,random_median_fx,num_essential,num_nonessential=resampleOneGo(go,go_terms,essential_go_counts,non_essential_go_counts,fx_df,n=n)
            rv_dict[go]=rv/n,this_median_fx,random_median_fx,num_essential,num_nonessential
    

    #adjust for multiple hypothesis

    rv_df=pd.DataFrame.from_dict(rv_dict, orient='index',columns=['raw_rv','go_median','median_in_random_sample','num_essential','num_nonessential'])

    pvalue_list=rv_df['raw_rv'].tolist()
    fdrbh_output = smm.multipletests(pvalue_list, method='fdr_bh') # benjamini hochberg method
    adjusted_pvalues = np.asarray(fdrbh_output[1].tolist())
		
    rv_df['bh_rv']=adjusted_pvalues
    
    rv_df.to_csv('resample_all_go_by_absolute_fxsize_n'+str(n)+'.tsv', sep='\t', index=True)

  
    return
    

if __name__ == "__main__":

    np.random.seed(777)
    random.seed(777)

    #parse external data files
    go_terms=ParseGO(go_file)
    essential_genes=ParseEssential(essential_file)

    #parse effect size from barseq
    fx_df=ParseFxSize(fxsize_file)

##    #shorten go_terms for debugging
##    go_terms=go_terms.head(n=100)

    print('...filtering...')
    

    #remove genes not in barseq from go terms
    go_terms=go_terms[go_terms['gene'].isin(fx_df['gene'].to_list())]

    #count
    

    all_go_counts=go_terms.groupby('go_term').count()
    all_go_counts.drop(columns=['gene_desc','gene'],inplace=True)
    all_go_counts.rename(columns={'sgd_name':'count'},inplace=True)
    #all_go_counts.to_csv('all_go_counts.tsv',sep=sep)
    print('there are '+str(len(all_go_counts))+' go terms in the barseq dataset')

    #filter go terms for ones with an acceptable number of genes

    
    acceptable_counts=all_go_counts[min_go_count<=all_go_counts['count']] #discard terms with too few genes
    print('there are '+str(len(acceptable_counts))+' go terms remaining after filtering out low gene counts(<'+str(min_go_count)+')')
    acceptable_counts=acceptable_counts[acceptable_counts['count']<=max_go_count]#discard terms wth too many genes
    print('there are '+str(len(acceptable_counts))+' go terms remaining after filtering out high gene counts(>'+str(max_go_count)+')')
    acceptable_count_go=acceptable_counts.reset_index()['go_term'].to_list()
    go_terms=go_terms[go_terms['go_term'].isin(acceptable_count_go)]


##    #filter for only bioprocess go terms #commented out to improve runtime because none of these changed
##    url='http://current.geneontology.org/ontology/go-basic.obo'
##    graph=obonet.read_obo(url)
##    bio_process=[]
##    for go in graph.nodes:
##        if graph.nodes[go]['namespace'] == 'biological_process':
##            bio_process += [go]
##    go_terms=go_terms[go_terms['gene'].isin(bio_process)]
##    print('there are '+str(len(acceptable_counts))+' go terms remaining after filtering for only bioprocess terms(>'+str(max_go_count)+')')

 
    

    #split go terms into essential and nonessential
    go_terms['essential']=go_terms['gene'].isin(essential_genes)
    essential_go_df=go_terms[go_terms['essential']==True]
    non_essential_go_df=go_terms[go_terms['essential']==False] 
    
    #count # essential and nonessential genes in each go term

    essential_go_counts=essential_go_df.groupby('go_term').count()
    essential_go_counts.drop(columns=['gene_desc','gene','essential'],inplace=True)
    essential_go_counts.rename(columns={'sgd_name':'count'},inplace=True)
   # essential_go_counts.reset_index(inplace=True)

    non_essential_go_counts=non_essential_go_df.groupby('go_term').count()
    non_essential_go_counts.drop(columns=['gene_desc','gene','essential'],inplace=True)
    non_essential_go_counts.rename(columns={'sgd_name':'count'},inplace=True)
   # non_essential_go_counts.reset_index(inplace=True)



 

  
    

    #resample
    resampleAllGo(go_terms,fx_df,essential_go_counts,non_essential_go_counts,n=n,low_n_file=low_n_file)
    

    
