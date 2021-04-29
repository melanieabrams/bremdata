import statsmodels.stats.multitest as smm
import sys
import pandas as pd
import numpy as np

LRT_in_spreadsheet='J3'

base_path='BSM_results/'

goi=['YGR198W', 'YMR207C', 'YGL082W', 'YNL049C', 'YDL035C', 'YDR508C', 'YBR136W',
 'YML099C', 'YPL254W', 'YIL152W', 'YKL017C', 'YGR140W', 'YJR127C', 'YDR375C',
 'YOR091W', 'YLR397C', 'YNL132W', 'YMR078C', 'YLR422W', 'YMR125W', 'YOR371C',
 'YMR094W', 'YMR167W', 'YDR103W', 'YDR318W', 'YAL026C', 'YDR180W', 'YOR092W',
 'YDR235W', 'YER151C', 'YMR275C', 'YKL114C', 'YOL081W', 'YPR049C', 'YGL095C',
 'YDR456W', 'YKL197C', 'YIL068C', 'YOR326W', 'YNR045W', 'YJR107W', 'YPL268W',
 'YJL062W', 'YCR042C'] #hits 37_2.0_1.1_..._10.0_3.0 with sc defect and effect size >0.5

#goi=['YGR198W']

#parse LRTs


#before this, put all bsm xls in a folder via these commands in Faisal's PAML directory:
##for folder in Y*; do echo $folder; cp ./$folder/BSM_* ./bsm_results/$folder$_BSM.xls; done #for top
##for folder in Y*; do echo $folder; cp ./$folder/results* ../bsm_results/$folder$_BSM.xls; done #for the rest
#MISSING: YBR136W
#for a couple of them, there were two files named BSM which threw an error, and I copied/named those manually


lrt_dict={}
missing_goi=[]
for gene in goi:

    

    filepath='r"'+gene+'.xls"'
    filepath=base_path+gene+'.xls'

    try:
        xls = pd.ExcelFile(filepath) #use r before absolute file path
        sheetX = xls.parse(0) #https://stackoverflow.com/questions/2942889/reading-parsing-excel-xls-files-with-python/50815107
        var1 = sheetX['Unnamed: 9']
        LRT=float(var1[1]) #this is the LRT p-value location in the xls output of Faisal's PAML
        lrt_dict[gene]=LRT
    except FileNotFoundError:
        missing_goi.append(gene)

    
    
print('missing:'+str(missing_goi))

print(lrt_dict)


df=pd.DataFrame.from_dict(lrt_dict, orient='index',columns=['lrt'])
#df=df.reset_index()

print(df)


###bh correction

pvalue_list=df['lrt'].tolist()
fdrbh_output = smm.multipletests(pvalue_list, method='fdr_bh') # benjamini hochberg method
adjusted_pvalues = np.asarray(fdrbh_output[1].tolist())
        
df['bh_lrt']=adjusted_pvalues

df.to_csv('bh_corrected_paml.tsv', sep='\t', index=True)
