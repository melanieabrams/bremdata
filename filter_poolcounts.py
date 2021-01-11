import pandas as pd
import sys

###Parameters

inFile=sys.argv[1]
#inFile='h1000_poolCount.txt'

cv_threshold=3.5

#conditions=['28C_for_36C','28C_for_37C','28C_for_38C','36C','37C','38C']#hardcoded in the col
condition_pairs={'36C':'28C_for_36C','37C':'28C_for_37C','38C':'28_for_38C'}

poolCounts = pd.read_csv(inFile, sep='\t') #parse infile


##RUN
#outfilename:
outFile=inFile.split('.')[0]+'_maxInsCV_'+str(cv_threshold)+'.txt'
logFile=outFile.split('.txt')[0]+'.log'

log_dict={}

#calculate cv for each condition

for condition in condition_pairs.keys():
    ref=condition_pairs[condition]
    
    condition_columns=[]
    ref_columns=[]
    
    condition_cv=condition+'_cv'
    ref_cv=ref+'_cv'
    
    for col in poolCounts.columns.values:
        if col.startswith(condition):
            condition_columns.append(col)
        if col.startswith(ref):
            ref_columns.append(col)
    poolCounts[condition_cv] = (poolCounts[condition_columns].std(axis=1)) / (poolCounts[condition_columns].mean(axis=1))
    poolCounts[ref_cv] = (poolCounts[condition_columns].std(axis=1)) / (poolCounts[condition_columns].mean(axis=1))

    #with both ref and condition
    #print(poolCounts[condition+'_cv'].mean(skipna=True)) #find mean cv: 2.0 to 2.7 to choose val
    #print(poolCounts[condition+'_cv'].std(skipna=True)) #find std of cv: 0.9 to 1.0
    #print(poolCounts[condition+'_cv'].max(skipna=True)) #max: 3.46
    #print(poolCounts[condition+'_cv'].min(skipna=True)) #min: 0.29
    
    #replace counts with 0s in both pairs if the insert has too high cv
    masked_from_condition = 0
    masked_from_ref = 0
    for index, row in poolCounts.iterrows(): #ex/ TATGGATCGTTGAGAGGGTC or GGAAGAACTACCCTTGTAAT at 38C (spot checked)
        if row[condition_cv]>cv_threshold:
            masked_from_condition+=1
            for column in condition_columns:
                row[column]=0     
            for column in ref_columns:
                row[column]=0 
            poolCounts.loc[index]=row
        elif row[ref_cv]>cv_threshold:
            masked_from_ref+=1
            for column in condition_columns:
                row[column]=0     
            for column in ref_columns:
                row[column]=0 
            poolCounts.loc[index]=row

    log_dict[condition]=masked_from_condition #log the number masked 
    log_dict[ref]=masked_from_ref

    #drop cv columns and save
    cvcols=[]
    for col in poolCounts.columns.values:
        if col.endswith('_cv'):
            cvcols.append(col)
    poolCounts=poolCounts.drop(columns=cvcols)
    poolCounts.to_csv(outFile,sep='\t',index=False)


#log fraction masked at ea condition
total_bc = len(poolCounts.index)
with open(logFile,'w') as wf:
    wf.writelines('filtered file: '+inFile+'\n')
    wf.writelines('maxInsCV: '+str(cv_threshold)+'\n')
    wf.writelines('condition pairs'+'\n')
    for condition in condition_pairs:
        ref=condition_pairs[condition]
        masked_condition=log_dict[condition]
        masked_ref=log_dict[ref]
        masked_c_pct=str(round(100.0*float((masked_condition)/float(total_bc)),3))+'%'
        masked_r_pct=str(round(100.0*float((masked_ref)/float(total_bc)),3))+'%'
        wf.writelines(condition+', with ref '+ref+': '+
                      masked_c_pct+' masked from condition, '+masked_r_pct+' masked from ref \n')

