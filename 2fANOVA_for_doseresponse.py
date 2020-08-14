from sys import argv
from os import path

import math
import numpy as np
import matplotlib.pyplot as plt

import csv

from scipy import stats
import pandas as pd
from patsy import dmatrices

import statsmodels.api as sm
from statsmodels.formula.api import ols
import statsmodels.stats.multitest as smm



##PARAMETERS##
infile="data_table_for_reorganization.csv"
outfile="esp1_dose_response_forRscript_unsorted.txt"
control_genotype="ScWT"


def parse_data():
    '''returns a list with [[genotype, temp, od] for all observations,
    and a list of all non-control genotypes'''
    
    all_data=[]
    test_gts=[]

    f=open(infile)
    next(f)
    for line in f:
        row_data=line.strip().split(',')
        #print(row_data)
        if row_data[0]=="Date":
            temps=row_data[3:]
            
        else:
            genotype=row_data[1]
            if genotype not in test_gts and genotype!=control_genotype:
                test_gts.append(genotype)
            biorep=row_data[1]+'_'+row_data[0]+'_'+row_data[2]
            ods=row_data[3:]
            for i in range(len(ods)):
                all_data.append([genotype,temps[i],float(ods[i])])
        
    
    f.close()
    return all_data, test_gts

def get_test_data(gt,all_data):
    '''returns the data that corresponds to the test strain + ctr'''
    test_data=[]
    for obs in all_data:
        if obs[0]==gt:
            test_data.append(obs)
        elif obs[0]==control_genotype:
            test_data.append(obs)
    return test_data
    

##START##
print("done loading libraries")
all_data,test_gts=parse_data()
print("done parsing data")

for gt in test_gts:
    print('==='+gt+'===')
    data=get_test_data(gt,all_data)
    
    data_df = pd.DataFrame(data, columns = ['STRAIN', 'TEMP','MEASUREMENT'])
    #print(data_df)

    species_temp_lm = ols('MEASUREMENT ~ C(STRAIN)*C(TEMP)',data_df).fit()
    anova_test = sm.stats.anova_lm(species_temp_lm,typ=2)
    pvalue=anova_test['PR(>F)'][2]
    
    Fval=anova_test['F'][2]
    print('Fvalue: '+str(Fval))
    print('Pvalue: '+str(pvalue))

