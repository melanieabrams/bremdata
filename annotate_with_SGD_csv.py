from sys import argv
import sys
import pandas as pd
import numpy as np
import math
import re

###USAGE###
#python3 annotate_with_SGD_csv.py sgdflatfile output_folder file_to_annotate

#SGD flatfile location: /Users/Melanie/Desktop/Melanie/NGS/SGD_features.tab

def parse_file(filename, sep='\t'): 
	df = pd.read_csv(filename, sep=sep)
	file_string = str(filename)
	p = re.compile('[^/]+[^\.](?=\.)') # matches everything before the period and after the last slash, ie. the identifier of the file
	m = p.search(file_string) # finds regex in file name
	file_identifier = m.group() # prints match in string format
	return df, file_identifier


def ParseSGDFlat(sgdflatfile):
    '''
    Parses SGD features flat file
    Input: SGD_features.tab file
    Output: 
    '''

    descr_dict={}
    
    f = open(sgdflatfile)
    for line in f:
        row_data =line.split("\t")
        yName=row_data[3]
        descr=row_data[15]
        descr_dict[yName]=descr
    return descr_dict



if __name__ == '__main__':
    
        sgdflatfile=sys.argv[1]
        output_folder = argv[2].strip('/')
        files = argv[3:]


        sgd_dict=ParseSGDFlat(sgdflatfile)l
        
        file_id_dict = {}
        df_list = []

        for each_file in files:
                df, file_identifier = parse_file(each_file)
                df.set_index('gene', inplace=True, drop=False)
                df_list.append(df)
                file_id_dict[file_identifier] = df

        for file_identifier, each_df in file_id_dict.items():
                gene_names=each_df['gene']
                annotations=[]
                for gene_name in gene_names:
                    annotations.append(sgd_dict[gene_name])
                each_df['Annotations'] = annotations
                each_df.to_csv(str(output_folder)+str('/')+str(file_identifier)+'.annotated_csv', sep='\t', index=False)

