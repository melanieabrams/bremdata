import seaborn as sns                                                
#paper_rc = {'lines.linewidth': 10}                  
#sns.set_context("paper", rc = paper_rc,font_scale=5)
#cmap = sns.cm.rocket_r
cmap='vlag'


import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

##USAGE: cumulative of the goi from 37 in 36C data

#TODO: get the other effect sizes. sort by sc_mean effect size?

##PARAMETERS
sc_mean_cutoff=-0.5

#fxsize files:
fx37='37C_2.0_1.1_coeffvarANDreadcutoff_uncollapsedMBA_unpaired_ratiocvfilt_1.6_3.0.all_fxsizes' #37C goi filtered from stats at the end of the pipeline
fx36='36C_2.0_1.1_coeffvarANDreadcutoff_uncollapsedMBA_unpaired_ratiocvfilt_1.6_3.0.all_fxsizes' #filtered 36

fxsizes={fx36:'36C',fx37:'37C'}
goi_from=fx37

#goi: 2_1.1_1.6_3.0, pcutoff 0.05, abs(sp_mean_obs)<=0.75, sc_mean_obs<-0.5
goi=['YGR198W','YLR397C','YBR136W','YGR140W','YNL049C','YKL017C','YGL082W',
     'YMR125W','YDR508C','YMR207C','YDR375C','YDR180W','YKL197C','YDR318W',
     'YMR094W','YOR326W','YBR081C','YPR049C','YIL152W','YER151C','YJR107W',
     'YAL026C','YDR456W','YLR141W','YPL268W','YDR235W']

def parse_fxsize(fxsize_filename, sep='\t'):
	'''parses the filtered gene inserts as a csv'''
	df = pd.read_csv(fxsize_filename,sep=sep,dtype={'effect_size':float,'sc_mean':float,'sp_mean':float})
	df = df.set_index('gene')
	return df

def filter_fxsize(fx_df,sc_mean_cutoff=sc_mean_cutoff):
	'''filters for goi'''
	print(fx_df.columns)
	filtered_df=fx_df[fx_df.index.isin(goi)]
	return filtered_df

def fxifingoi(gene,fx):
	'''helper for lambda'''
	if gene=='YMR207C':
		print('HI')
		print(fx)
	if gene in goi:
		return fx
	else:
		return np.nan

def transform_data(fx_df,fxsize):
	'''drops unneeded columns and calculates nonabsolute fx size.  adds a column which only has fxsize if in goi'''
	#drop all but new 'one sided fx' size and gene name
	fx_df['nonabsolute_fx']=fx_df['sp_mean']-fx_df['sc_mean']
	all_col=fx_df.columns.tolist()
	col_to_keep=['nonabsolute_fx','gene']
	col_to_drop=[col for col in all_col if col not in col_to_keep]
	dropped=fx_df.drop(columns=col_to_drop)
	dropped=dropped.rename(columns={'nonabsolute_fx':fxsize})

	dropped=dropped.reset_index()


	dropped['36Cdata_37Cgoi']=dropped.apply(lambda x: fxifingoi(x['gene'],x[fxsize]),axis=1)
	print(dropped)

	return dropped


def plotCDF(fx_df):
	#print(fx_df)
	fig, ax = plt.subplots(figsize=(6,20))
	mask = fx_df.isnull()
	ax=sns.kdeplot(data=fx_df.filter(like="36C", axis="columns"),common_norm=False,multiple='stack') #for histograms
	#ax=sns.kdeplot(data=fx_df.filter(like="36C", axis="columns"),cumulative=True, common_norm=False, common_grid=True) #for estimated cdf
	ax.set_xlabel('Effect Size')
	plt.tick_params(axis='both', which='major', labelsize=10, labelbottom = True, bottom=False, top = False, labeltop=False)

##	median_36=fx_df['36C'].mean()
##	median_36data_37Cgoi=fx_df['36Cdata_37Cgoi'].mean()
##
##	plt.vlines(median_36,0,1.5)
##	plt.vlines(median_36data_37Cgoi,0,1.5)
	
	plt.show()
	#plt.savefig(figName)

	
def getCDFInput(fxsize):
	condition=fxsizes[fxsize]
	fxsize_df=parse_fxsize(fxsize)
	#filtered_data=filter_fxsize(fxsize_df,sc_mean_cutoff=sc_mean_cutoff)
	#transformed_data=transform_data(filtered_data,condition)
	transformed_data=transform_data(fxsize_df,condition)
	return transformed_data

###RUN###

if __name__ == "__main__":
	fxsize_HeatMapInput=[]
	for fxsize in fxsizes:
		if goi_from!=fxsize:
			plotCDF(getCDFInput(fxsize))

			


