from sys import argv
import sys
import pandas as pd
import numpy as np
import math
import re

### USAGE ###
# python fitness_ratios.py output_folder file_location/filtered_file.filtered_inserts
# for each insert, calculates log2(39/28) using final normalized read counts


def parse_file(filename, sep='\t'): 
        df = pd.read_csv(filename, sep=sep)
        file_string = str(filename)
        p = re.compile('[^/]+[^\.](?=\.)') # matches everything before the period and after the last slash, ie. the identifier of the file
        m = p.search(file_string) # finds regex in file name
        file_identifier = m.group() # prints match in string format
        return df, file_identifier


if __name__ == '__main__':

        output_folder = argv[1].strip('/')
        files = argv[2:]

        file_id_dict = {}
        df_list = []

        for each_file in files:
                df, file_identifier = parse_file(each_file)
                df.set_index('ID', inplace=True, drop=False)
                df_list.append(df)
                file_id_dict[file_identifier] = df

        for file_identifier, each_df in file_id_dict.iteritems():
                columns = each_df.columns.values
                read_columns = [col for col in columns if '_averaged_reads'or '_n_av' in col]
                reads_28 = []
                reads_39 = []
                reads_T0 = []
                for each_col in read_columns:
                        if each_col.startswith('28'):
                                reads_28.append(each_col)
                        if each_col.startswith('39'):
                                reads_39.append(each_col)
                        if each_col.startswith('T0'):
                                reads_T0 = each_col
                
                #add median column
                brs_28= [br for br in reads_28 if '_n_av' in br]
                median28=each_df.groupby(['28_br_averaged_reads'])[brs_28].apply(np.nanmedian)
                print(median28)
                median28.name='MEDIAN'
                median_df=each_df.join(median28, on=['28_br_averaged_reads'])
                
                for biorep39 in reads_39:
                        repID = biorep39[2:]
                        
                        div39_28 = median_df[biorep39]/median_df['MEDIAN']
                        colName = '39_28_log2' + repID
                        each_df[colName] = np.log2(div39_28)
                                                               
                #read_columns = [col for col in columns if '_averaged_reads' in col] CHANGED FROM THIS TO ALLOW PRESERVATION BIOREP STR
##              for each_col in read_columns:
##                      if each_col.startswith('28'):
##                              reads_28 = each_col
##                      if each_col.startswith('39'):
##                              reads_39 = each_col
##                      if each_col.startswith('T0'):
##                              reads_T0 = each_col
##              div39_28 = each_df[reads_39] / each_df[reads_28]
##              each_df['39_28_log2'] = np.log2(div39_28)
                each_df.to_csv(str(output_folder)+str('/')+str(file_identifier)+'_unpaired_byMedian.insert_ratios', sep='\t', index=False)
