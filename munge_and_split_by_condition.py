import pandas as pd
#Parameters

inFileName='geneFit_all_1readfilter_individual_replicates.txt'
#outFileName='geneFit_all_for_Tstat_MA.txt'
outFileName='default'
conditions=['36C','37C','38C']
#note: set to work with unpaired data, so form of naming convention is 38C_11_vs_38C_averaged_references for biorep 11 at 38C

#Functions


def mungeAndSplit(inFileName,conditions,outFileName=outFileName):
    '''input: geneFit_individual_replicates.txt file from RBSeq
    output: geneFit_temp_for_Tstat file.txt file for downstream analysis - one for all, one split by condition'''
    
    #make one munged file, save
    df = pd.read_csv(inFileName, sep='\t') #parse infile
    df = df[df.NearestGene != 'intergenic'] #remove intergenic locations
    df = df.drop(labels=['NearestGene','Annotation'],axis=1) #drop columns not used in _for_Tstat outfile
    df = df.rename(columns={'AlternateID':'annotation'}) #keep annotation in form scYHR023W, rename col
    df['gene']=df['annotation'].str.slice(start=2)
    df['allele']=df['annotation'].str.slice(stop=2)
    if outFileName=='default':
        outFileName=inFileName.split('_individual_replicates')[0]+'_forTstat_MA.txt' #save outfile not split by temp
    df.to_csv(outFileName, sep='\t', index=False)

    #split by condition and save
    for condition in conditions:
        columns_to_drop = []
        renamed_columns = {}
        for column in df.columns:
            if 'averaged_references' in column:
                if condition not in column:
                    columns_to_drop.append(column)
                else:
                    br_n=column.split('_')[1]
                    br_T=column.split('_')[0]
                    br_ref=column.split('_vs')[1]
                    renamed_columns[column]='Biorep_'+br_n+'_'+br_T+'_vs'+br_ref
        condition_df=df.drop(columns=columns_to_drop)
        condition_df=condition_df.rename(columns=renamed_columns)
        conditionOutFileName=condition+'_'+outFileName
        print(condition_df)
        condition_df.to_csv(conditionOutFileName, sep='\t', index=False)

    return None

##RUN##

mungeAndSplit(inFileName,conditions, outFileName=outFileName)
                    
                
                

