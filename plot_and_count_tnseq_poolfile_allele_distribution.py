import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

###PARAMETERS###
poolfile='/opt/shared/skerker/rhseq/rhseq_barseq/ReferenceFiles/piggyBAC_all20_poolfile_sep1_carly.txt'

###FUNCTIONS###
def parse_poolfile(poolfile):
    '''parse poolfile df and add allele column'''
    poolfile_df=pd.read_csv(poolfile,sep='\t')
    poolfile_df['allele']=poolfile_df.apply(lambda x: str(x['scaffold'])[:2],axis=1)
    #poolfile_df=poolfile_df.head(10000) # for testing, b/c plotting long time
    return poolfile_df



def count_alleles(poolfile_df,allele):
    '''count alleles matching poolfile'''
    df=poolfile_df[poolfile_df['allele']==allele]
    bc_count=len(df.barcode.to_list())
    return bc_count

###RUN###

poolfile_df=parse_poolfile(poolfile)
print(poolfile_df.head())
print(poolfile_df.columns)

sc_alleles=count_alleles(poolfile_df,'sc')
print('sc: '+str(sc_alleles))

sp_alleles=count_alleles(poolfile_df,'sp')
print('sp: '+str(sp_alleles))

#plot
print('plotting')

sc_abundances=poolfile_df[poolfile_df['allele']=='sc']['n'].to_list()
sp_abundances=poolfile_df[poolfile_df['allele']=='sp']['n'].to_list()
fig, ax = plt.subplots(figsize=(9,6))
print('here')
sns.histplot(x=sc_abundances, label='Insert in cerevisiae',color='#F1B629')
sns.histplot(x=sp_abundances, label='Insert in paradoxus',color='#08A5CD')
print('here')

plt.legend()
ax.set_xlabel('n')

plt.savefig('poolfile_readcounts_hist.png')
ax.set_xscale('log')
plt.savefig('poolfile_readcounts_hist_logscale.png')

