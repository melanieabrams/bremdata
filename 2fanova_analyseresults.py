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

#to use python3.7:  #source activate /Applications/anaconda3

##PARAMETERs

#data as list of lists: [strain, temp, fraction sm+med cells]
data=[[62,28,0.122807018],[62,28,0.222222222],
      [68,28,0.240384615],[68,28,0.197368421],
      [98,28,0.282051282],[98,28,0.257575758],
      [62,39,0.054794521],[62,39,0.111111111],
      [68,39,0.147368421],[68,39,0.152173913],
      [98,39,0.022727273],[98,39,0.094827586]]

#data list of lists: [strain, temp, fraction sm+med cells] for only 68 and 98
data=[[68,28,0.240384615],[68,28,0.197368421],
      [98,28,0.282051282],[98,28,0.257575758],
      [68,39,0.147368421],[68,39,0.152173913],
      [98,39,0.022727273],[98,39,0.094827586]]




####  Running the script  ####



data_df = pd.DataFrame(data, columns = ['STRAIN', 'TEMP','MEASUREMENT'])

species_temp_lm = ols('MEASUREMENT ~ C(STRAIN)*C(TEMP)',data_df).fit()
anova_test = sm.stats.anova_lm(species_temp_lm,typ=2)
pvalue=anova_test['PR(>F)'][2] 
Fval=anova_test['F'][2]
print('Fvalue: '+str(Fval))
print('Pvalue: '+str(pvalue))
