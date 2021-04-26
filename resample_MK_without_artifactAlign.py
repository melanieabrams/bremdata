import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns



import sys
sys.path.append('/usr2/people/mabrams/scripts/')
from resample_cdubin import resample_med


#modified from Claire Dubin's script, found at https://github.com/clairedubin/thermotolerance/blob/master/MK_analysis.ipynb

##PARAMETERS###
gff='/usr2/people/mabrams/Amended_Genomes/S288C/saccharomyces_cerevisiae_R64-1-1_20110208_withoutFasta.gff'

figname='MK_barcoded_rhseq_hits.png'

essential_file='/usr2/people/clairedubin/sacc/external_datasets/essential.csv' #essential genes - from Winzeler et al., 1999

goi=['YGR198W','YLR397C','YBR136W','YGR140W','YKL017C','YGL082W',
     'YMR125W','YDR508C','YMR207C','YDR375C','YDR180W','YKL197C','YDR318W',
     'YMR094W','YOR326W','YBR081C','YPR049C','YIL152W','YER151C','YJR107W',
     'YDR456W','YLR141W','YPL268W','YDR235W']
    #without YNL049C (SFB2) and YAL026C(DRS2) which had artifactual alignments not beginning with start codon 

##FUNCTIONS###

def ParseFromGFF(gff):
    '''
    Parses gff
    Output: dict of {yName:{geneName}}
    '''
    ann_dict={}
    
    f = open(gff)
    lines=[]
    for line in f:
        if line[0]!='#': #skip header rows
            row_data=line.split("\t")
            if row_data[2]=='gene':
                info=row_data[8].split(";")
                yName=info[0].split('=')[1]
                geneName=info[2].split('=')[1]
                ann_dict[yName]=geneName
    f.close()
    return ann_dict


###RUN###

#parse essential file
essential = pd.read_csv(essential_file, header=None)
essential[1] = essential[1].str.strip('\t')
essential_genes = [i.split(' ')[0] for i in essential[1]]


#build gene dict
ann_dict=ParseFromGFF(gff)

gene_dict = {}
for yName in goi:
    geneName=ann_dict[yName]
    gene_dict[yName]=geneName

#read MK input, calculate NI
df = pd.read_csv('/usr2/people/clairedubin/sacc/current_sacc_datasets/MK_1011Scer_EuroSparNASparB_031920.csv')
df.head()

df['Dn/Ds'] = df['Dn']/df['Ds']
df['Pn/Ps'] = df['Pn']/df['Ps']
df['NI'] = df['Pn/Ps']/df['Dn/Ds']

df['alpha'] = 1 - (df['Ds']*df['Pn'])/(df['Dn']*df['Pn'])


df['NI'].median()


candidates = df[df['gene'].isin(gene_dict.keys())]
candidates['name'] = [gene_dict[g] for g in candidates['gene']]
candidates[['gene','name', 'Ds', 'Dn', 'Ps', 'Pn', 'Dn/Ds', 'Pn/Ps', 'NI', 'alpha',]]

##Resample NI
np.random.seed(7777)

print('missing: ', [gene_dict[i] for i in gene_dict.keys() if i not in candidates['gene'].tolist()])
print('NI p={}'.format(resample_med(candidates, df, 'NI', essential_genes, direction='less_than', graph=True,n=10000,figName='resampling_distribution_MK.png')))

#plot NI by gene

plt.figure(figsize=(10, 5))
plt.bar([gene_dict[gene] for gene in candidates.gene], candidates['NI'])
plt.ylabel("NI")
plt.tick_params(axis='x', labelbottom=True, labeltop=False, bottom=False)
plt.xticks(rotation=90)
plt.axhline(y=candidates['NI'].median(),linewidth=1, color='orange', label='candidate median')
plt.axhline(y=df['NI'].median(),linewidth=1, color='purple', label='genome wide median')

plt.legend(loc='upper right', bbox_to_anchor=(1, 1.1),
      ncol=3, fancybox=True, shadow=True)
plt.savefig(figname)
                                                    

                                                      
