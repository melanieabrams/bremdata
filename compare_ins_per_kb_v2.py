import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import random

###PARAMETERS###


fgi='1.5_1.1_coeffvarANDreadcutoff_uncollapsedMBA_unpaired_ratiocvfilt_20.0_2.0.filtered_gene_inserts'
mwu='1.5_1.1_coeffvarANDreadcutoff_uncollapsedMBA_unpaired_ratiocvfilt_20.0_2.0.mwu_test_results'
gff='saccharomyces_cerevisiae_R64-1-1_20110208_withoutFasta.gff'
save_intermediate=True #to save a csv of the usable ins per gene
##STATUS: resampling not implemented, but saved usable ins and kb count df###

dxy_euro_file='dxy_1011pops_EuroSpar.csv'

goi=['YLR397C','YGR098C','YMR168C','YKR054C',
     'YHR023W','YDR180W','YPL174C','YCR042C',
     'YMR016C','YJR135C','YJL025W','YDR443C',
     'YKL134C'] #from popgen MS

essential_file='essential.csv' #source: http://www-sequence.stanford.edu/group/yeast_deletion_project/Essential_ORFs.txt


##FUNCTIONS###

def ParseFromGFF(gfffile):
    '''
    Parses SGD features flat file
    Input: SGD_features.tab file
    Output: dict of gene and gene length
    '''
    
    gene_length={}
    
    f = open(gfffile)
    for line in f:
        if line[0]!='#' : #skip header rows
            row_data=line.split("\t")
            annotation_type=row_data[2]
            if annotation_type=='CDS':
                chrom=row_data[0]
                start=int(row_data[3])
                stop=int(row_data[4])
                length=abs(start-stop)
                info=row_data[8].split(";")
                yName=info[0].split('=')[1]
                if yName[0]=='Y' and len(yName)>5:
                    gene_length[yName]=length
                
    f.close()
    
    return gene_length

def ParseEssential(essential_file):
    '''load list of essential genes as essential_genes'''
    print('...parsing essential genes...')
    essential = pd.read_csv(essential_file, header=None)
    essential[1] = essential[1].str.strip('\t')
    essential_genes = essential[1].str.strip(' ').tolist()
    return essential_genes

def parseFGI(filename,sep='\t'):
    '''return the count for each gene'''
    df = pd.read_csv(filename,sep=sep)
    #print(df[df['gene']=='YPL174C'])

    print('parsed '+filename)
    return df

def parseMWU(mwu,gene_length_dict):
    '''get dict of gene:avg g1'''
    p_dict={}
    f = open(mwu)
    next(f) #skip header
    for line in f:
        row_data=line.split("\t")
        p_dict[row_data[0]]=float(row_data[2])
    for gene in gene_length_dict: #avoid key errors by passing nan for missing data
        if gene not in p_dict:
            p_dict[gene]=np.nan
    return p_dict

def parseDxyEuro(dxy_file,gene_length_dict):
    dxy_euro= pd.read_csv(dxy_euro_file, header=None) #sp euro
    dxy_euro.columns = ['population', 'gene', 'dxy', 'spar_strain_count', 'scer_strain_count']
    dxy_euro = dxy_euro[dxy_euro['spar_strain_count'] >= 8]
    #print(dxy_euro.head)
    dxy_wine=dxy_euro[dxy_euro['population']=='1WineEuropean'] #sc euro too
    #print(dxy_wine)
    dxy_dict=pd.Series(dxy_wine.dxy.values,index=dxy_wine.gene).to_dict()
    for gene in gene_length_dict: #avoid key errors by passing nan for missing data
        if gene not in dxy_dict:
            dxy_dict[gene]=np.nan
    return dxy_dict
    

def calcInsPerKb(df,dxy_dict,p_dict,gene_length_dict):
    df=df['gene'].value_counts().rename_axis('gene').to_frame('usable_ins').reset_index()#calc usable ins per gene
    df=df[df['usable_ins']>0]
    print(df)
    df['length_bp']=df.apply(lambda x:gene_length_dict[x['gene']],axis=1)
    df['usable_ins_per_kb']=df.apply(lambda x:(float(x['usable_ins'])/float(x['length_bp'])),axis=1)
    df=df[df['length_bp']>100]
    print(df)
    df['dxy']=df.apply(lambda x:dxy_dict[x['gene']],axis=1)
    df['adj_p']=df.apply(lambda x:p_dict[x['gene']],axis=1)
    df.to_csv('usable_ins_and_length_and_dxy_and_mwu.csv')

    print(df)
    return df


def plotScatter(df):
    #all x and y
    x=df['usable_ins_per_kb'].to_list()
    y=df['dxy'].to_list()

    #just goi x and y (for red dots)
    df['ingoi']=df.apply(lambda x:x['gene'] in goi,axis=1)
    goi_df=df[df['ingoi']==True]
    goi_x=goi_df['usable_ins_per_kb'].to_list()
    goi_y=goi_df['dxy'].to_list()


            
    line1=plt.scatter(x,y, color='black')
    line2=plt.scatter(goi_x,goi_y, color='red')
    #plt.legend((line1, line2), ('sp', 'sc'))

    plt.ylabel('dxy (1WineEuro sc vs. Euro sp)')
    plt.xlabel('usable_ins_per_kb')
    plt.rc('xtick',labelsize=80)
    plt.rc('ytick',labelsize=80)
    plt.show()



def get_goi_df(df):
    print('...filtering for goi')

    #filter for goi
    df['ingoi']=df.apply(lambda x: x['gene'] in goi,axis=1)
    df=df[df.ingoi==True]
    df=df.drop(columns=['ingoi'])
    print('goi in 37C dataset: '+str(len(df)) +'out of '+str(len(goi)))
    print(df['gene'].to_list())
    print()
    #df = df.set_index('gene')

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

        fgi_df=parseFGI(fgi)
        gene_length_dict=ParseFromGFF(gff)

        p_dict=parseMWU(mwu,gene_length_dict)
        dxy_dict=parseDxyEuro(dxy_euro_file,gene_length_dict)
        
        df=calcInsPerKb(fgi_df,dxy_dict,p_dict,gene_length_dict)
        plotScatter(df)

        
##        goi_df=get_goi_df(df)
##    
##
##       
##        #print(all37)
##
##        essential_genes=ParseEssential(essential_file)
##
##    
##        #split goi into essential and non essential genes
##
##        essential_goi = [gene for gene in goi if gene in essential_genes]
##        essential_count=len(essential_goi)
##        print('\nessential goi ('+str(essential_count)+'):')
##        print(essential_goi)
##        non_essential_goi=[gene for gene in goi if gene not in essential_genes]
##        non_essential_count=len(non_essential_goi)
##        print('\nnonessential goi ('+str(non_essential_count)+'):')
##        print(non_essential_goi)
##
##        #resample
##        resample(df36,goi36,n=10000)

