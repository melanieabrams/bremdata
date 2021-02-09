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
fx37='2.0_1.1_coeffvarANDreadcutoff_uncollapsedMBA_unpaired_ratiocvfilt_1.6_3.0pcut0pt05spmaxoffby0pt75.fxsizes' #37C goi filtered from stats at the end of the pipeline
fx36='36C_2.0_1.1_coeffvarANDreadcutoff_uncollapsedMBA_unpaired.all_fxsizes' #36 and 38C all fx sizes from step 7 logratio files, look at goi within these
fx38='38C_2.0_1.1_coeffvarANDreadcutoff_uncollapsedMBA_unpaired.all_fxsizes'

sortby=fx38
fxsizes={fx36:'36C',fx37:'37C',fx38:'38C'}

gen_36=[12.84584039,12.96574443,12.82881066,12.47694361,12.69612127,12.57667477,12.72712914,12.42038924,12.4493805,12.52735329,13.08141686,13.31136996]
gen_37=[10.71381697,10.6059797,10.70591248,10.74969728,10.72307974,10.87017294,10.95210186,10.61683005,9.636340693,10.18646952,10.72492607,10.55877579]
gen_38=[5.123086751,5.10433666,5.121015401,5.203592714,5.037821465,5.022367813,4.83541884,4.882643049,4.832890014,4.794415866,4.860466259,4.87036472]

avg_gen={'36C':np.mean(gen_36),'37C':np.mean(gen_37),'38C':np.mean(gen_38)} #right now, just norm to the average - will need to change script 10effect_size if I wanna do by each br


def parse_fxsize(fxsize_filename, sep='\t'):
	'''parses the filtered gene inserts as a csv'''
	df = pd.read_csv(fxsize_filename,sep=sep,dtype={'effect_size':float,'sc_mean':float,'sp_mean':float})
	df = df.set_index('gene')
	##	if sortby==fxsize_filename:
##			df.sort_values(by=['sc_mean'],ascending=True,inplace=True)
	return df

def filter_fxsize(fx_df,sc_mean_cutoff=sc_mean_cutoff):
	'''filters for goi'''
	filtered_df=fx_df[fx_df['sc_mean']<-0.5]
	return filtered_df


def transform_and_norm_data(fx_df,fxsize):
	'''drops unneeded columns, normalizes to gen'''
	#drop all but fx size and gene name
	all_col=fx_df.columns.tolist()
	col_to_keep=['sp_mean','gene']
	col_to_drop=[col for col in all_col if col not in col_to_keep]
	dropped=fx_df.drop(columns=col_to_drop)
	dropped=dropped.rename(columns={'sp_mean':fxsize})
	dropped[fxsize]=dropped[fxsize].div(avg_gen[fxsize])
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
	transformed_data=transform_and_norm_data(filtered_data,condition)
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

			


