import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

###PARAMETERS###
poolfile='piggyBAC_all20_poolfile_sep1_carly.txt'
poolcounts='poolCount_37C.txt' # the 37C one: renamed from poolCount.txt

###FUNCTIONS###
def parse_poolfile(poolfile):
    '''parse poolfile df and add allele column and reverse bc and rcbarcode to match poolcounts'''
    poolfile_df=pd.read_csv(poolfile,sep='\t')
    poolfile_df['allele']=poolfile_df.apply(lambda x: str(x['scaffold'])[:2],axis=1)
    poolfile_df.drop(columns=['barcode'],inplace=True)
    poolfile_df.rename(columns={'rcbarcode':'barcode'},inplace=True)
    #poolfile_df=poolfile_df.head(10000) # for testing, b/c plotting long time
    return poolfile_df

def parse_poolcounts(poolcounts):
    '''parse poolfile df and add allele column and a combined counts column'''
    poolcounts_df=pd.read_csv(poolcounts,sep='\t')
    poolcounts_df['allele']=poolcounts_df.apply(lambda x: str(x['scaffold'])[:2],axis=1)
    temp_col=[col for col in poolcounts_df.columns if 'C_' in col]
    poolcounts_df['poolcount_n']=poolcounts_df[temp_col].sum(axis=1)
   # poolcounts_df=poolcounts_df.head(10000) # for testing, b/c plotting long time)
    return poolcounts_df

def name_xypair(n,poolcount_n):
    #print(str(n)+','+str(poolcount_n))
    return str(n)+','+str(poolcount_n)

def get_merged_df(poolfile_df,poolcounts_df):
    '''merge, and count pairs of counts for each barcode'''


    merged_df=pd.merge(poolfile_df,poolcounts_df, on='barcode')
    print(merged_df)
    #merged_df['xy']=merged_df['n']+merged_df['poolcount_n']
    merged_df['xy']=merged_df.apply(lambda x: name_xypair(x['n'],x['poolcount_n']),axis=1)
    merged_df['xycount']=merged_df.groupby('xy')['xy'].transform('count')
    print(merged_df)
    print(merged_df.columns)
    return merged_df

###RUN###

poolfile_df=parse_poolfile(poolfile)
poolcounts_df=parse_poolcounts(poolcounts)
print(poolfile_df.head())
print(poolcounts_df.head())
merged_df=get_merged_df(poolfile_df,poolcounts_df)

fig=plt.figure()
x=merged_df['n']
y=merged_df['poolcount_n']
s=merged_df['xycount']

plt.scatter(x=x,y=y,s=s,facecolors='none',edgecolors='black')
plt.xscale('log')
plt.yscale('log')
plt.xlabel('Tn-Seq counts')
plt.ylabel('BarSeq counts (across replicates)')
plt.xlim(1,2*max(x))
plt.ylim(1,2*max(y))
plt.show()

