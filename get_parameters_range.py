import pandas as pd
import numpy as np
import collections

pd.options.mode.chained_assignment = None  # default='warn'

##Parameters###
scorefile='blinded_scores_dyad_size_parsable.csv'
strain_key={'62_39C':['A','C'],'68_39C':['E','G'],'98_39C':['I','K'], # key to blinded images
            '62_28C':['B','D'],'68_28C':['F','H'],'98_28C':['J','L']}
outfilename='blinded_scores_dyad_size_parsable_params.txt'


###Functions###


def ParseBlindedScores(scorefile,strain_key,sep=','):
    '''load blinded scores, add a column for all dyads scored, strain, and the image set'''
    df = pd.read_csv(scorefile,sep=sep,dtype={'lg':int,'med':int,'sm':int}) #read the dataframe
    df['sum']=df['lg']+df['med']+df['sm'] # get the sum of dyads counted in each image
    df['image_set']=df.apply(lambda x:str(x['File']).split('-')[0],axis=1) #get the image set it came from

    image_set_dict={}
    for strain in strain_key:
        for image_set in strain_key[strain]:
            image_set_dict[image_set]=strain
    df['strain']=df.apply(lambda x:image_set_dict[x['image_set']],axis=1)
    #print(df)
    return df

def calcParams(df):
    '''calculate range and median for image deck, write outfile'''
    strains=set(df['strain'].to_list())
    strain_params={} #key: strain. values: [min, max, median,mean,stdev,total]
    for strain in strains:
        strain_df=df[df['strain']==strain]
        #print(strain_df)
        strain_sum=strain_df['sum'].to_list()
        strain_min=min(strain_sum)
        strain_max=max(strain_sum)
        strain_median=np.median(strain_sum)
        strain_mean=np.mean(strain_sum)
        strain_std=np.std(strain_sum)
        strain_total=sum(strain_sum)
        strain_params[strain]=[strain_min,strain_max,strain_median,strain_mean,strain_std,strain_total]
    print(strain_params)

    with open(outfilename,'w') as wf:
        wf.writelines('strain\tmin\tmax\tmedian\tmean\tstdev\ttotal\n')
        for strain in strains:
            wf.writelines(strain+'\t'+'\t'.join([str(i) for i in strain_params[strain]])+'\n')
        

    return
    


    
    


if __name__ == "__main__":

    #parse external data files
    scores=ParseBlindedScores(scorefile,strain_key)
    calcParams(scores)
    
    
