import pandas as pd
import numpy as np
from scipy import stats
import statsmodels.stats.multitest as smm
from funcy import flatten, isa

###PARAMETERS###
obs36='36C_2.0_1.1_coeffvarANDreadcutoff_uncollapsedMBA_unpaired.all_fxsizes'
obs37='37C_2.0_1.1_coeffvarANDreadcutoff_uncollapsedMBA_unpaired_ratiocvfilt_1.6_3.0.all_fxsizes'

save_goi_obs=True #set to False if you don't want to write intermediate files of just the goi observations
savename36='goi_'+obs36
savename37='goi_'+obs37

#goi: 2_1.1_1.6_3.0, pcutoff 0.05, abs(sp_mean_obs)<=0.75, sc_mean_obs<-0.5
goi=['YGR198W','YLR397C','YBR136W','YGR140W','YNL049C','YKL017C','YGL082W',
     'YMR125W','YDR508C','YMR207C','YDR375C','YDR180W','YKL197C','YDR318W',
     'YMR094W','YOR326W','YBR081C','YPR049C','YIL152W','YER151C','YJR107W',
     'YAL026C','YDR456W','YLR141W','YPL268W','YDR235W']

###FUNCTIONS###
def parseInsRatios(filename,sep='\t'):
    df = pd.read_csv(filename,sep=sep,dtype={'effect_size':float})
    df=df.drop(columns=['sc_mean','sp_mean','mwu_pval'])
    print('parsed '+filename)
    return df

def get_goi_logr(df,savename,save=True):
    print('...filtering for goi')

    #filter for goi
    df['ingoi']=df.apply(lambda x: x['gene'] in goi,axis=1)
    df=df[df.ingoi==True]
    df=df.drop(columns=['ingoi'])
    df = df.set_index('gene')

    #save, print, and return
    print(df)
    if save==True:
        df.to_csv(savename, sep='\t', index=True)
    return df

def getAllVals(df):
    valLists=df.values.tolist()
    flatList=list(flatten(valLists))
    cleanList=[i for i in flatList if str(i)!='nan']
    #print(cleanList[:50])
    return cleanList

def mwu(list1,list2):
    mwu_test = stats.mannwhitneyu(list1,list2)   # mwu test
    two_tailed_pval = mwu_test[1] * 2.0
    return two_tailed_pval

###RUN###
if __name__ == "__main__":
        logr36=parseInsRatios(obs36)
        goi_logr_36=get_goi_logr(logr36,savename36,save=save_goi_obs)
        
        logr37=parseInsRatios(obs37)          
        goi_logr_37=get_goi_logr(logr37,savename37,save=save_goi_obs)

        all36=getAllVals(goi_logr_36)
        all37=getAllVals(goi_logr_37)

        mwu_36and37=mwu(all36,all37)

        print(mwu_36and37) #result: 1.292590369227476e-07
