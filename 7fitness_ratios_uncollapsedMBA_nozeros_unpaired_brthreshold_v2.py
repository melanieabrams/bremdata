from sys import argv
import sys
import pandas as pd
import numpy as np
import math
import re

### USAGE ###
# python fitness_ratios.py output_folder file_location/filtered_file.filtered_inserts
# for each biorep of each insert above a threshold, calculates log2(39/28) using final normalized read counts

#PARAMETERS
br_read_thresh = 35.0 #threshold for number of reads for a 39C biorep or 28C average must have before it has a logratio calculated

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
                #get columns corresponding to bioreps or averages across bioreps 
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

        
                #replace zeros with nans so not taken into the average
                each_df.replace(to_replace=0.0, value=np.nan, inplace=True) 


                
                for biorep39 in reads_39:
                        repID = biorep39[2:]
                        for biorep in reads_28:
                                if biorep == '28_br_averaged_reads':
                                        avg28 = biorep
                        aboveThresholdBR39=each_df[biorep39].where(each_df[biorep39]>br_read_thresh,np.nan)
                        aboveThresholdAvg28=each_df[avg28].where(each_df[avg28]>br_read_thresh,np.nan)
                        #print(aboveThreshold39)
                        div39_28 = aboveThresholdBR39/aboveThresholdAvg28
                        #div39_28 = each_df[biorep39]/each_df[avg28]
                        colName = '39_28_log2' + repID
                        print(colName)
                        each_df[colName] = np.log2(div39_28)



                #each_df=each_df.replace([np.inf, -np.inf,np.nan], 'None')
                each_df=each_df.replace([np.inf, -np.inf], np.nan)
                print('writing csv...')
		each_df.to_csv(str(output_folder)+str('/')+str(file_identifier)+'_unpaired_'+str(int(br_read_thresh))+'avgandbrrt.insert_ratios', sep='\t', index=False)
