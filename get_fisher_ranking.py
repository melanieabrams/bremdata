from sys import argv
import sys
import pandas as pd
import numpy as np
import math
import re

###USAGE###
#python3 get_fisher_ranking.py fisher_test output_folder file_to_annotate

#SGD flatfile location: /Users/Melanie/Desktop/Melanie/NGS/SGD_features.tab

def parse_file(filename, sep='\t'): 
	df = pd.read_csv(filename, sep=sep)
	file_string = str(filename)
	p = re.compile('[^/]+[^\.](?=\.)') # matches everything before the period and after the last slash, ie. the identifier of the file
	m = p.search(file_string) # finds regex in file name
	file_identifier = m.group() # prints match in string format
	return df, file_identifier


def rankFisher(fisherfile):
    descr_dict={}

    lines_counter=1
    
    f = open(fisherfile)
    for line in f:
        row_data =line.split("\t")
        yName=row_data[0]
        descr_dict[yName]=lines_counter
        lines_counter+=1
    return descr_dict



if __name__ == '__main__':
    
        fisherfile=sys.argv[1]
        output_folder = argv[2].strip('/')
        files = argv[3:]


        fisher_dict=rankFisher(fisherfile)
        
        file_id_dict = {}
        df_list = []

        for each_file in files:
                df, file_identifier = parse_file(each_file)
                df.set_index('gene', inplace=True, drop=False)
                df_list.append(df)
                file_id_dict[file_identifier] = df

        for file_identifier, each_df in file_id_dict.items():
                gene_names=each_df['gene']
                fisher_rank=[]
                for gene_name in gene_names:
                        try:
                                fisher_rank.append(fisher_dict[gene_name])
                        except KeyError:
                                fisher_rank.append('unranked')
                each_df['fisher_rank'] = fisher_rank
                each_df.to_csv(str(output_folder)+str('/')+str(file_identifier)+'.ranked_csv', sep='\t', index=False)

