import seaborn as sns                                                
#paper_rc = {'lines.linewidth': 10}                  
#sns.set_context("paper", rc = paper_rc,font_scale=5)
#cmap = sns.cm.rocket_r

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

##USAGE: make a simple heatmap (one effect size per condition)

#TODO: get the other effect sizes. sort by sc_mean effect size?

##PARAMETERS
sc_mean_cutoff=-0.5

#fxsize files:
fx37='37C_2.0_1.1_coeffvarANDreadcutoff_uncollapsedMBA_unpaired_ratiocvfilt_1.6_3.0.all_fxsizes' #37C goi filtered from stats at the end of the pipeline
fx36='36C_2.0_1.1_coeffvarANDreadcutoff_uncollapsedMBA_unpaired.all_fxsizes' #36 and 38C all fx sizes from step 7 logratio files, look at goi within these
fx38='38C_2.0_1.1_coeffvarANDreadcutoff_uncollapsedMBA_unpaired.all_fxsizes'

sortby=fx38
fxsizes={fx36:'36C',fx37:'37C',fx38:'38C'}

goi=['YLR397C','YGR098C','YMR168C','YKR054C',
     'YHR023W','YDR180W','YPL174C','YCR042C',
     'YMR016C','YJR135C','YJL025W','YDR443C',
     'YKL134C'] #from popgen MS

def parse_fxsize(fxsize_filename, sep='\t'):
	'''parses the filtered gene inserts as a csv'''
	df = pd.read_csv(fxsize_filename,sep=sep,dtype={'effect_size':float,'sc_mean':float,'sp_mean':float})
	df = df.set_index('gene')
	##	if sortby==fxsize_filename:
##			df.sort_values(by=['sc_mean'],ascending=True,inplace=True)
	return df

def filter_fxsize(fx_df,sc_mean_cutoff=sc_mean_cutoff):
	'''filters for goi'''
	print(fx_df.columns)
	filtered_df=fx_df[fx_df.index.isin(goi)]
	return filtered_df


def transform_data(fx_df,fxsize):
	'''drops unneeded columns, normalizes to gen'''
	#drop all but fx size and gene name
	all_col=fx_df.columns.tolist()
	col_to_keep=['sc_mean','gene']
	col_to_drop=[col for col in all_col if col not in col_to_keep]
	dropped=fx_df.drop(columns=col_to_drop)
	dropped=dropped.rename(columns={'sc_mean':fxsize})
	#dropped['filler']=1
	#print(dropped)

	#pivot
	#pivoted=dropped.pivot(index='None',columns=['effect_size'],values='effect_size')

	#return pivoted
	return dropped

def plotHeatmap(fx_df):
	p=sns.heatmap(fx_df)
	plt.show()

	
def getHeatMapInput(fxsize):
	condition=fxsizes[fxsize]
	fxsize_df=parse_fxsize(fxsize)
	filtered_data=filter_fxsize(fxsize_df,sc_mean_cutoff=sc_mean_cutoff)
	transformed_data=transform_data(filtered_data,condition)
	#plotHeatmap(transformed_data)
	return transformed_data

###RUN###

if __name__ == "__main__":
	fxsize_HeatMapInput=[]
	for fxsize in fxsizes:
		fxsize_HeatMapInput.append(getHeatMapInput(fxsize))

	mergedHeatMapInput=pd.concat(fxsize_HeatMapInput, axis=1, join="inner")
	mergedHeatMapInput.sort_values(by=[fxsizes[sortby]],ascending=True,inplace=True)
	print(mergedHeatMapInput)
	plotHeatmap(mergedHeatMapInput)

			

