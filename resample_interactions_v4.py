import pandas as pd
import numpy as np
import collections
import random
import statsmodels.stats.multitest as smm

pd.options.mode.chained_assignment = None  # default='warn'

##Parameters###

intx_type='any' #default: 'any'
intx_file='interaction_data.tab' #source: http://sgd-archive.yeastgenome.org/curation/literature/interaction_data.tab ; downloaded 2/19/2021
essential_file='essential.csv' #source: http://www-sequence.stanford.edu/group/yeast_deletion_project/Essential_ORFs.txt

goi=['YGR198W','YLR397C','YBR136W','YGR140W','YNL049C','YKL017C','YGL082W',
     'YMR125W','YDR508C','YMR207C','YDR375C','YDR180W','YKL197C','YDR318W',
     'YMR094W','YOR326W','YBR081C','YPR049C','YIL152W','YER151C','YJR107W',
     'YAL026C','YDR456W','YLR141W','YPL268W','YDR235W']

carlygoi=['YLR397C','YGR098C','YKR054C',
     'YHR023W','YDR180W','YPL174C','YCR042C',
     'YMR016C','YJR135C','YJL025W','YDR443C',
     'YKL134C'] #from popgen MS, excluding YMR168C b/c only one br

goi=goi+carlygoi

#this one analyzes the number of intxns between over number of total intxns

#note: this go term enrichment script is modified from Dubin et al., 2020 (https://github.com/clairedubin/peroxisome_evolution/blob/68589c20c2e637d645ab4e79511db3ea992b121f/gene_expression/expression_analysis_artieri.ipynb)

###Functions###

##def CompliantYNames(bait,hit):
##    if CompliantYName(bait)==True:
##        if CompliantYName(hit)==True:
##            return True
def CompliantYName(sgdname):
    if sgdname.startswith('Y'):
         if sgdname.endswith('W') or sgdname.endswith('C'):
            return True

def ParseInteractions(intx_file,sep='\t'):
    '''load interaction data as intx_terms'''
    print('...parsing interaction file...')
    intx_terms = pd.read_csv(intx_file, header=None,sep=sep,dtype=str)
    intx_terms = intx_terms.drop(columns=[1,3,4,6,7,8,9,10,11])
  
    intx_terms = intx_terms.rename(columns={0: 'bait', 2:'hit', 5:'intx_type'})
 
    intx_terms['bait'] = [i[0] for i in intx_terms['bait'].str.split('|')]
    intx_terms['hit'] = [i[0] for i in intx_terms['hit'].str.split('|')]

    intx_terms = intx_terms.drop_duplicates()


    intx_terms['compliant_bait']=intx_terms.apply(lambda x:CompliantYName(x['bait'])==True,axis=1)
    intx_terms['compliant_hit']=intx_terms.apply(lambda x:CompliantYName(x['hit'])==True,axis=1)

    intx_terms=intx_terms[(intx_terms['compliant_hit']==True)&(intx_terms['compliant_hit']==True)]
    intx_terms=intx_terms.drop(columns=['compliant_bait','compliant_hit'])
    
    
    #print(intx_terms[intx_terms['bait']=='YPR200C'])
    return intx_terms

def ParseEssential(essential_file):
    '''load list of essential genes as essential_genes'''
    print('...parsing essential genes...')
    essential = pd.read_csv(essential_file, header=None)
    essential[1] = essential[1].str.strip('\t')
    essential_genes = essential[1].str.strip(' ').tolist()
    return essential_genes



def resampleAllInteractions(goi,intx_terms,essential_count,non_essential_count,n=10000,save_intermediate=True,intx_type='any'):
    if intx_type!='any':
        intx_terms=intx_terms[intx_terms['intx_type']==intx_type]
        intx_terms=intx_terms.drop_duplicates()
        
    print('...resampling all interactions...')

    #drop all self-interactions 
    intx_terms['self-intxn']=intx_terms.apply(lambda x:x['bait']==x['hit'],axis=1)
    intx_terms=intx_terms[intx_terms['self-intxn']==False]
    intx_terms=intx_terms.drop(columns=['self-intxn'])
    ##TRY HERE###


    
    #split interactions into essential and nonessential bait, and if bait in goi

    intx_terms['essential']= intx_terms['bait'].isin(essential_genes)
    intx_terms['bait_in_goi']=intx_terms ['bait'].isin(goi)
    
    essential_intx_df=intx_terms[intx_terms['essential']==True]
    non_essential_intx_df=intx_terms[intx_terms['essential']==False]
    
    goi_intx=intx_terms[intx_terms['bait_in_goi']==True]

    print('- saving intermediate?: '+str(save_intermediate))
    if save_intermediate==True:
        goi_intx.to_csv('goi_intx_df.tsv', sep='\t', index=True)
    

    #count interactions where both bait and hit are in goi
    goi_any_count=len(goi_intx.bait.tolist())
    print('- number of interactions involving any if the goi: '+str(goi_any_count))
    goi_goi_intx=goi_intx[goi_intx.hit.isin(goi)] #hit in goi too
    goi_goi_count=len(goi_goi_intx)
    print('- number of interactions where the bait is also in the goi:'+str(goi_goi_count))
    goi_goi_ratio=float(goi_goi_count)/float(goi_any_count)
    print('- ratio of (goi_goi) to (goi_any):'+str(goi_goi_ratio))
    print(goi_goi_intx)
    

    ##get n random samples of the same # of essential and nonessential genes
    print('- number of random samples: '+str(n))
    random_samples=[]
    sample_dfs=[]
    for i in range(n):
        random_sample = list(random.sample(essential_intx_df.bait.tolist(),essential_count))
        random_sample+= list(random.sample(non_essential_intx_df.bait.tolist(),non_essential_count))
        #random_samples.append(random_sample)
        intx_terms['bait_in_random']=intx_terms['bait'].isin(random_sample)
        intx_terms['hit_in_random']=intx_terms['hit'].isin(random_sample)
        sample_any_df=intx_terms[intx_terms['bait_in_random']==True]
        sample_sample_df=intx_terms[(intx_terms['bait_in_random']==True)&(intx_terms['hit_in_random']==True)]
        sample_dfs.append(sample_sample_df)
        sample_sample_ratio=float(len(sample_sample_df))/float(len(sample_any_df))
        random_samples.append(sample_sample_ratio)
        intx_terms=intx_terms.drop(columns=['bait_in_random','hit_in_random'])
    #print(sample_dfs[-1])
    print(random_samples)
    

    rv=0.0
    for samp_samp_ratio in random_samples:
        if samp_samp_ratio>=goi_goi_ratio:
            rv+=1
    
    print('resampling value: '+str(rv/n))
            
  
    return
    

if __name__ == "__main__":

    #parse external data files
    intx_terms=ParseInteractions(intx_file)
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

    print('interaction type: '+intx_type)



    #resample
    resampleAllInteractions(goi,intx_terms,essential_count,non_essential_count,n=10000,save_intermediate=False,intx_type=intx_type)
    

    
