from sys import argv
import sys
import pandas as pd
import numpy as np
import math
import re

### USAGE ###
# python fitness_ratios.py output_folder file_location/filtered_file.filtered_inserts
# for each insert, calculates log2(exptl/ctrl) using final normalized read counts


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

        for file_identifier, each_df in file_id_dict.items():
                each_df.replace(to_replace=0.0, value=np.nan, inplace=True)
                columns = each_df.columns.values
                read_columns = [col for col in columns if '_averaged_reads'or '_n_av' in col]
                reads_ctrl = []
                reads_exptl = []
                reads_T0 = []
                for each_col in read_columns:
                        if each_col.startswith('ctrl'):
                                reads_ctrl.append(each_col)
                        if each_col.startswith('exptl'):
                                reads_exptl.append(each_col)
                        if each_col.startswith('T0'):
                                reads_T0 = each_col
                for biorepexptl in reads_exptl:
                        repID = biorepexptl[2:]
                        for biorep in reads_ctrl:
                                if biorep == 'ctrl_br_averaged_reads':
                                        avgctrl = biorep
                        divexptl_ctrl = each_df[biorepexptl]/each_df[avgctrl]
                        colName = 'exptl_ctrl_log2' + repID
                        print(colName)
                        each_df[colName] = np.log2(divexptl_ctrl)



                #replaces ratios where either numerator or denominator is zero with nan 
                each_df=each_df.replace([np.inf, -np.inf], np.nan)
                print('writing csv...')
                each_df.to_csv(str(output_folder)+str('/')+str(file_identifier)+'_unpaired.insert_ratios', sep='\t', index=False)
