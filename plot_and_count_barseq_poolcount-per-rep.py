import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

###PARAMETERS###
poolcounts='poolCount_37C.txt' # the 37C one: renamed from poolCount.txt

###FUNCTIONS###
def parse_poolcounts(poolcounts):
    '''parse poolfile df and add allele column and a column of combined counts'''
    poolcounts_df=pd.read_csv(poolcounts,sep='\t')
    poolcounts_df['allele']=poolcounts_df.apply(lambda x: str(x['scaffold'])[:2],axis=1)
    temp_col=[col for col in poolcounts_df.columns if 'C_' in col]
    poolcounts_df['poolcount_n']=poolcounts_df[temp_col].sum(axis=1)
   # poolcounts_df=poolcounts_df.head(10000) # for testing, b/c plotting long time)
    return poolcounts_df


def count_alleles(poolcounts_df,allele):
    '''count alleles matching poolcounts'''
    df=poolcounts_df[poolcounts_df['allele']==allele]
    bc_count=len(df.barcode.to_list())
    return bc_count

def count_biorep(poolcounts_df,rep):
    '''return counts of bc present in one biorep'''
    rep_list=poolcounts_df[rep].to_list() #get counts for that rep
    clean_list = [x for x in rep_list if np.isnan(x) == False and x!=0] #remove nan or 0
    clean_count=len(clean_list) # get len
    return clean_count 

def count_bioreps(poolcount_df):
    '''return counts of bc present in all bioreps'''
    temp_col=[col for col in poolcounts_df.columns if 'C_' in col]
    temp_col_counts=[]
    for col in temp_col:
        br_count=count_biorep(poolcount_df,col)
        print(col+': '+str(br_count))
        temp_col_counts.append(br_count)
    return temp_col_counts
        
        

###RUN###

#parse
poolcounts_df=parse_poolcounts(poolcounts)

#count alleles
sc_alleles=count_alleles(poolcounts_df,'sc')
print('sc: '+str(sc_alleles))

sp_alleles=count_alleles(poolcounts_df,'sp')
print('sp: '+str(sp_alleles))

#count bioreps for each allele
sc_df=poolcounts_df[poolcounts_df['allele']=='sc']
sp_df=poolcounts_df[poolcounts_df['allele']=='sp']
print('counting sc')
sc_counts=count_bioreps(sc_df)
print('counting sp')
sp_counts=count_bioreps(sp_df)
print('sc_counts,sp_counts by biorep:')
print(sc_counts,sp_counts)

#plot
fig, ax = plt.subplots(figsize=(9,6))
print('here')
sns.histplot(x=sc_counts, label='Insert in cerevisiae',color='#F1B629')
sns.histplot(x=sp_counts, label='Insert in paradoxus',color='#08A5CD')
print('here')

plt.legend()
ax.set_xlabel('BC in a replicate tube')
plt.savefig('poolcounts_replicate_readcounts_dist.png')

#plot overall allele distribution


print('plotting all')

sc_abundances=poolcounts_df[poolcounts_df['allele']=='sc']['poolcount_n'].to_list()
sp_abundances=poolcounts_df[poolcounts_df['allele']=='sp']['poolcount_n'].to_list()

fig, ax = plt.subplots(figsize=(9,6))
print('here1')
sns.histplot(x=sc_abundances, label='Insert in cerevisiae',color='#F1B629')
sns.histplot(x=sp_abundances, label='Insert in paradoxus',color='#08A5CD')
print('here1')

plt.legend()
ax.set_xlabel('n')

plt.savefig('poolcounts_readcounts_hist.png')
ax.set_xscale('log')
plt.savefig('poolcounts_readcounts_hist_logscale.png')

