import pandas as pd
from sys import argv
import sys
import numpy as np

import matplotlib.pyplot as plt
import seaborn as sns

### USAGE ###
# python plot_hist_cv.py

##PARAMS##

files=['28.normalized_averaged_bioreps','39.normalized_averaged_bioreps']

coeff_var_cutoff = 2.0 #cutoff for coefficient of variation
norm_reads_cutoff = 1.1 # cutoff for normalized read count
filter_strategy = 'both' # can be 'coeff' or 'reads' or 'both'



def parse_file(filename, sep='\t'): 
        df = pd.read_csv(filename, sep=sep)
        return df

global test
test = 0.0

def fillna_and_combine(grp):
        # keeps track of how many inserts have gone through, prints out so you can see how long the script is taking to run
        global test
        test += 1
        if test % 10000 == 0:
                print(test)

        grp['ID'] = grp.index
        filled = grp.fillna(method='ffill').fillna(method='bfill')
        return filled.iloc[0]

def make_and_populate_new_columns(grp):
        
        condition = grp['condition'].iloc[0] #first row of 0th column
        #condition = grp.ix[0, 'condition'] #deprecated, previous
    

        if condition != 'T0':
                column_name = condition+'_br_averaged_reads'
                grp[column_name] = grp['br_averaged_reads']
                column_name2 = condition+'_br_coeffvar'
                grp[column_name2] = grp['br_coeffvar']
        return grp


if __name__ == '__main__':

        


        file_dict = {}
        df_list = []

        for each_file in files:
                df = parse_file(each_file)
                df = df.head(1000) #THIS LINE SHOULD BE COMMENTED OUT AFTER TEST RUN
                df.set_index('ID', inplace=True, drop=False)
                condition_ID = df['condition'][1]
                if condition_ID == 'T0':
                        df.drop('bio_rep_ID', axis=1, inplace=True)
                file_dict[condition_ID] = df
                df_list.append(df)
        
        concatd_df = pd.concat(df_list, axis=0)
        all_columns = list(concatd_df.columns.values)
        drop_columns = [col for col in all_columns if '_norm_' in col]
        drop_columns.extend([col for col in all_columns if '_tr_' in col])
        concatd_df.drop(drop_columns, axis=1, inplace=True)
        concatd_df.condition = concatd_df.condition.astype(str)
        concatd_df.set_index('condition', inplace=True, drop=False)
        concatd_df.index.name = None # required for later versions of python than 2.6
        grouped2 = concatd_df.groupby('condition', sort=False).apply(make_and_populate_new_columns)
        grouped2.set_index('ID', inplace=True, drop=False)
        grouped2.drop(['br_averaged_reads', 'br_coeffvar', 'condition'], axis=1, inplace=True)
        grouped2['NEWID'] = grouped2['ID']
        grouped_tmp = grouped2.groupby('NEWID', sort=False)
        group3 = grouped_tmp.apply(fillna_and_combine)

        group3.rename(columns={'39_br_coeffvar':'37_br_coeffvar',"39_br_averaged_reads":"37_br_averaged_reads"},inplace=True) #for plot
        columns = list(group3.columns.values)

        
#       check why a whole bunch None here
        group3['sameval']=group3.apply(lambda x: x['37_br_averaged_reads']==x['28_br_averaged_reads'],axis=1)
        group3[group3['sameval']==True].to_csv('sameval_diagnosing_step6.tsv',sep='\t')

        
        
        
        #columns_to_plot = [col for col in columns if 'br_coeffvar' in col] # filter by either of the "conditions" coeffvar
        columns_to_plot = [col for col in columns if 'br_averaged_reads' in col] # filter by either of the "reads" coeffvar
        for column in columns_to_plot:
            fig, ax = plt.subplots(figsize=(6,6))
            print(group3)
            print(group3[column])
            cv_rows=group3[~group3[column].eq('None')]
            print(cv_rows[column])
            sns.histplot(x=cv_rows[column].astype(float))
            #plt.axvline(coeff_var_cutoff,color='red')
            plt.axvline(norm_reads_cutoff,color='red',ls='--',lw=0.5)
            ax.set_xscale('log')
            plt.savefig(column+'_hist.png')

        
        

        
