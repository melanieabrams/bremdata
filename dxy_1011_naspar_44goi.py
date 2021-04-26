import pandas as pd
import numpy as np
import csv

#modified from https://github.com/clairedubin/thermotolerance/blob/master/dxy_analysis_1011pops_NASpar.ipynb


##PARAMETERS###

n=10000 #n=10k default

goi=['YGR198W', 'YMR207C', 'YGL082W', 'YNL049C', 'YDL035C', 'YDR508C', 'YBR136W',
 'YML099C', 'YPL254W', 'YIL152W', 'YKL017C', 'YGR140W', 'YJR127C', 'YDR375C',
 'YOR091W', 'YLR397C', 'YNL132W', 'YMR078C', 'YLR422W', 'YMR125W', 'YOR371C',
 'YMR094W', 'YMR167W', 'YDR103W', 'YDR318W', 'YAL026C', 'YDR180W', 'YOR092W',
 'YDR235W', 'YER151C', 'YMR275C', 'YKL114C', 'YOL081W', 'YPR049C', 'YGL095C',
 'YDR456W', 'YKL197C', 'YIL068C', 'YOR326W', 'YNR045W', 'YJR107W', 'YPL268W',
 'YJL062W', 'YCR042C'] #hits 37_2.0_1.1_..._10.0_3.0 with sc defect and effect size >0.5

#essentialfile='/usr2/people/clairedubin/publishable_data/external_datasets/essential.csv'
essentialfile='essential.csv'
dxyfile='dxy_1011pops_NASpar.csv'

outfile='dxy_1011_naspar_resample_37C_hits_44goi.csv'


pop_dict={'10FrenchGuianahuman': 31,
 '11Alebeer': 18,
 '12WestAfricancocoa': 13,
 '13Africanpalmwine': 28,
 '14CHNIII': 2,
 '15CHNII': 2,
 '16CHNI': 1,
 '17Taiwanese': 3,
 '18FarEastAsia': 9,
 '19Malaysian': 6,
 '1WineEuropean': 268,
 '1WineEuropeansubclade1': 18,
 '1WineEuropeansubclade2': 13,
 '1WineEuropeansubclade3': 24,
 '1WineEuropeansubclade4': 39,
 '20CHNV': 2,
 '21Ecuadorean': 10,
 '22FarEastRussian': 4,
 '23NorthAmericanoak': 13,
 '24Asianislands': 11,
 '25Sake': 47,
 '26Asianfermentation': 39,
 '2Alpechin': 17,
 '3Brazilianbioethanol': 35,
 '4Mediterraneanoak': 8,
 '5Frenchdairy': 32,
 '6Africanbeer': 20,
 '7Mosaicbeer': 21,
 '8Mixedorigin': 72,
 '9Mexicanagave': 7,
 'M1Mosaicregion1': 17,
 'M2Mosaicregion2': 20,
 'M3Mosaicregion3': 113}


##FUNCTIONS
def resample(df1, df2, n=10000, graph=False):

    
    results = 0

    actual_med = df1['dxy'].median()
    
    essential_count = len([i for i in df1['gene'] if i in essential_genes])
    nonessential_count = len(df1['gene']) - essential_count
    
    print('candidate gene median dxy: ', actual_med)
    print('essential count: '+str(essential_count)+', nonessential_count: '+str(nonessential_count))
    print('resampling pool size: ', df2.shape[0])
    
    essential_df = df2[df2['gene'].isin(essential_genes)]
    nonessential_df =  df2[~df2['gene'].isin(essential_genes)]

    sample_meds = []
    
    for _ in range(n):
        sample = essential_df.sample(n=essential_count)
        sample = sample.append(nonessential_df.sample(n=nonessential_count))
        sample_med = sample['dxy'].median()
        sample_meds += [sample_med]
        if sample_med >= actual_med:
            results += 1
    
    if graph:
        
        print('p=',results/n)
        sns.distplot(sample_meds)
        plt.axvline(x=actual_med, color='r', label = 'True median')
        plt.xlabel('Sample medians')
        plt.ylabel('Frequency')
    
    return results/n

###RUN###

if __name__ == "__main__":

    #essential genes + GO term info
    essential = pd.read_csv(essentialfile, header=None)
    essential[1] = essential[1].str.strip('\t')
    essential_genes = essential[1].str.strip(' ').tolist()


    #load raw dxy data
    all_dxy = pd.read_csv(dxyfile, header=None) ##CHANGE
    all_dxy.columns = ['population', 'gene', 'dxy', 'spar_strain_count', 'scer_strain_count']


    #drop any rows where spar_strain_count < 8 or scer_strain count < 75% of the population

    all_dxy = all_dxy[all_dxy['spar_strain_count'] >= 8]
    all_dxy.shape

    #resample

    print('n='+str(n))

    np.random.seed(777)

    p_dict = {}

    for pop in pop_dict:
        size = pop_dict[pop]
        
        df = all_dxy[all_dxy['population']==pop]
        df = df[df['scer_strain_count']>= .75*size]
        
        candidates = df[df['gene'].isin(goi)]
        
        print('')
        print('---- {} ----'.format(pop))

        print('missing: ', [i for i in goi if i not in candidates['gene'].tolist()])
        
        p = resample(candidates, df,n=n)
        p_dict[pop] = [df['dxy'].median(), candidates['dxy'].median(), df.shape[0], p]
        
        print('p = ', p)

    p_df = pd.DataFrame.from_dict(p_dict, orient="index",columns=['genome_median','goi median','genes tested','rv'])
    p_df.to_csv(outfile)
