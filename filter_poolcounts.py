import pandas as pd
import sys

###Parameters

#inFile=sys.argv[1]
inFile='h1000_poolCount.txt'
cv_threshold=3.4

conditions=['28C_for_36C','28C_for_37C','28C_for_38C','36C','37C','38C']#hardcoded in the col
condition_pairs={'36C':'28C_for_36C','37C':'28C_for_37C','38C':'28_for_38C'}

poolCounts = pd.read_csv(inFile, sep='\t') #parse infile


##RUN
#outfilename:
outFile=inFile.split('.')[0]+'_maxInsCV_'+str(cv_threshold)+'.txt'


#calculate cv for each condition

for condition in condition_pairs.keys():
    ref=condition_pairs[condition]
    condition_columns=[]
    ref_columns=[]
    for col in poolCounts.columns.values:
        if col.startswith(condition):
            condition_columns.append(col)
        if col.startswith(ref):
            ref_columns.append(col)
    poolCounts[condition+'_cv'] = (poolCounts[condition_columns].std(axis=1)) / (poolCounts[condition_columns].mean(axis=1))
    poolCounts[ref+'_cv'] = (poolCounts[condition_columns].std(axis=1)) / (poolCounts[condition_columns].mean(axis=1))

    #with both ref and condition
    #print(poolCounts[condition+'_cv'].mean(skipna=True)) #find mean cv: 2.0 to 2.7 to choose val
    #print(poolCounts[condition+'_cv'].std(skipna=True)) #find std of cv: 0.9 to 1.0
    #print(poolCounts[condition+'_cv'].max(skipna=True)) #max: 3.46
    #print(poolCounts[condition+'_cv'].min(skipna=True)) #min: 0.29
    
    #replace counts with 0s in both pairs if the insert has too high cv
    for index, row in poolCounts.iterrows(): #ex/ TATGGATCGTTGAGAGGGTC or GGAAGAACTACCCTTGTAAT at 38C (spot checked)
        if row[condition+'_cv']>cv_threshold:
            for column in condition_columns:
                row[column]=0     
            for column in ref_columns:
                row[column]=0 
            poolCounts.loc[index]=row

    #drop cv columns and save
    cvcols=[]
    for col in poolCounts.columns.values:
        if col.endswith('_cv'):
            cvcols.append(col)
    poolCounts=poolCounts.drop(columns=cvcols)
    poolCounts.to_csv(outFile,sep='\t',index=False)
