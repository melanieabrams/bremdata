import pandas as pd
import numpy as np
import random

###PARAMETERS###



obs36='36C_2.0_1.1_coeffvarANDreadcutoff_uncollapsedMBA_unpaired.all_fxsizes' #unfiltered 36C genes
obs36='36C_2.0_1.1_coeffvarANDreadcutoff_uncollapsedMBA_unpaired_ratiocvfilt_1.6_3.0.all_fxsizes' #old filter scheme 36C data

save_goi_obs=True #set to False if you don't want to write intermediate files of just the goi observations
savename36='goi_'+obs36

goi=['YGR198W', 'YMR207C', 'YGL082W', 'YNL049C', 'YDL035C', 'YDR508C', 'YBR136W',
 'YML099C', 'YPL254W', 'YIL152W', 'YKL017C', 'YGR140W', 'YJR127C', 'YDR375C',
 'YOR091W', 'YLR397C', 'YNL132W', 'YMR078C', 'YLR422W', 'YMR125W', 'YOR371C',
 'YMR094W', 'YMR167W', 'YDR103W', 'YDR318W', 'YAL026C', 'YDR180W', 'YOR092W',
 'YDR235W', 'YER151C', 'YMR275C', 'YKL114C', 'YOL081W', 'YPR049C', 'YGL095C',
 'YDR456W', 'YKL197C', 'YIL068C', 'YOR326W', 'YNR045W', 'YJR107W', 'YPL268W',
 'YJL062W', 'YCR042C'] #hits 37_2.0_1.1_..._10.0_3.0 with sc defect and effect size >0.5


essential_file='essential.csv' #source: http://www-sequence.stanford.edu/group/yeast_deletion_project/Essential_ORFs.txt


##FUNCTIONS###

def ParseEssential(essential_file):
    '''load list of essential genes as essential_genes'''
    print('...parsing essential genes...')
    essential = pd.read_csv(essential_file, header=None)
    essential[1] = essential[1].str.strip('\t')
    essential_genes = essential[1].str.strip(' ').tolist()
    return essential_genes

def parseInsRatios(filename,sep='\t'):
    df = pd.read_csv(filename,sep=sep,dtype={'sc_mean':float,'sp_mean':float})
    print(df[df['gene']=='YPL174C'])
        
    df['directional_effect_size']=df['sp_mean']-df['sc_mean']
    #df=df.drop(columns=['effect_size','mwu_pval','sc_mean','sp_mean'])
    print('parsed '+filename)
    return df

def get_goi_df(df,savename,save=True):
    print('...filtering for goi')

    #filter for goi
    df['ingoi']=df.apply(lambda x: x['gene'] in goi,axis=1)
    df=df[df.ingoi==True]
    df=df.drop(columns=['ingoi'])
    print('goi in 37C dataset: '+str(len(df)) +'out of '+str(len(goi)))
    print(df['gene'].to_list())
    print()
    #df = df.set_index('gene')

    #save, print, and return
    print(df)
    if save==True:
        df.to_csv(savename, sep='\t', index=True)
    return df


def resample(df,goi_df,n=10000):

    goi_median=goi_df['directional_effect_size'].median()

    df['essential']= df['gene'].isin(essential_genes)
    
    essential_df=df[df['essential']==True]
    non_essential_df=df[df['essential']==False]

    ##get n random samples of the same # of essential and nonessential genes
    print('- number of random samples: '+str(n))
    random_samples=[]
    for i in range(n):
        random_sample = list(random.sample(essential_df.gene.tolist(),essential_count))
        random_sample+= list(random.sample(non_essential_df.gene.tolist(),non_essential_count))
        random_df=df[df['gene'].isin(random_sample)==True]
        random_df_median=random_df['directional_effect_size'].median()
        random_samples.append(random_df_median)
    #print(sample_dfs[-1])
    print(random_samples)
    

    rv=0.0
    for random_median in random_samples:
        if random_median>=goi_median:
            rv+=1
    p=(rv/float(n))
    print(p)
    return p



###RUN###
if __name__ == "__main__":

        df36=parseInsRatios(obs36)          
        goi36=get_goi_df(df36,savename36,save=save_goi_obs)
    

       
        #print(all37)

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

        #resample
        resample(df36,goi36,n=10000)

