import seaborn as sns                                                
#paper_rc = {'lines.linewidth': 10}                  
#sns.set_context("paper", rc = paper_rc,font_scale=5)
#cmap = sns.cm.rocket_r
cmap='vlag'


import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from matplotlib import rcParams
rcParams.update({'figure.autolayout': True})

##USAGE: make a simple heatmap (one effect size per condition)

#TODO: get the other effect sizes. sort by sc_mean effect size?

##PARAMETERS
sc_mean_cutoff=-0.5

#fxsize files:
fx37='37C_2.0_1.1_coeffvarANDreadcutoff_uncollapsedMBA_unpaired_ratiocvfilt_10.0_3.0pcut0pt05spmaxoffby100.fxsizes' #37C goi filtered from stats at the end of the pipeline
fx36='36C_2.0_1.1_coeffvarANDreadcutoff_uncollapsedMBA_unpaired.all_fxsizes' #36 and 38C all fx sizes from step 7 logratio files, look at goi within these
#fx38='38C_2.0_1.1_coeffvarANDreadcutoff_uncollapsedMBA_unpaired.all_fxsizes'

#fx37='37C_2.0_1.1_coeffvarANDreadcutoff_uncollapsedMBA_unpaired.all_fxsizes' #unfiltered 37
#fx36='36C_2.0_1.1_coeffvarANDreadcutoff_uncollapsedMBA_unpaired_ratiocvfilt_1.6_3.0.all_fxsizes' #filtered 36

sortby=fx37
#fxsizes={fx36:'36C',fx37:'37C',fx38:'38C'}
fxsizes={fx36:'36C',fx37:'37C'}


gen_36=[12.84584039,12.96574443,12.82881066,12.47694361,12.69612127,12.57667477,12.72712914,12.42038924,12.4493805,12.52735329,13.08141686,13.31136996]
gen_37=[10.71381697,10.6059797,10.70591248,10.74969728,10.72307974,10.87017294,10.95210186,10.61683005,9.636340693,10.18646952,10.72492607,10.55877579]
gen_38=[5.123086751,5.10433666,5.121015401,5.203592714,5.037821465,5.022367813,4.83541884,4.882643049,4.832890014,4.794415866,4.860466259,4.87036472]

avg_gen={'36C':np.mean(gen_36),'37C':np.mean(gen_37),'38C':np.mean(gen_38)} #right now, just norm to the average - will need to change script 10effect_size if I wanna do by each br


goi=['YGR198W', 'YMR207C', 'YGL082W', 'YNL049C', 'YDL035C', 'YDR508C', 'YBR136W',
 'YML099C', 'YPL254W', 'YIL152W', 'YKL017C', 'YGR140W', 'YJR127C', 'YDR375C',
 'YOR091W', 'YLR397C', 'YNL132W', 'YMR078C', 'YLR422W', 'YMR125W', 'YOR371C',
 'YMR094W', 'YMR167W', 'YDR103W', 'YDR318W', 'YAL026C', 'YDR180W', 'YOR092W',
 'YDR235W', 'YER151C', 'YMR275C', 'YKL114C', 'YOL081W', 'YPR049C', 'YGL095C',
 'YDR456W', 'YKL197C', 'YIL068C', 'YOR326W', 'YNR045W', 'YJR107W', 'YPL268W',
 'YJL062W', 'YCR042C'] #hits 37_2.0_1.1_..._10.0_3.0 with sc defect and effect size >0.5


gff='saccharomyces_cerevisiae_R64-1-1_20110208_withoutFasta.gff'

###FUNCTIONS$$$
def ParseFromGFF(gff):
	'''
	Parses gff
	Output: dict of {yName:{geneName}}
	'''
	ann_dict={}
	
	f = open(gff)
	lines=[]
	for line in f:
		if line[0]!='#': #skip header rows
			row_data=line.split("\t")
			if row_data[2]=='gene':
				info=row_data[8].split(";")
				yName=info[0].split('=')[1]
				geneName=info[2].split('=')[1]
				ann_dict[yName]=geneName
	f.close()
	return ann_dict



def parse_fxsize(fxsize_filename, sep='\t'):
	'''parses the filtered gene inserts as a csv'''
	df = pd.read_csv(fxsize_filename,sep=sep,dtype={'effect_size':float,'sc_mean':float,'sp_mean':float})
	df = df.set_index('gene')
	##      if sortby==fxsize_filename:
##                      df.sort_values(by=['sc_mean'],ascending=True,inplace=True)
	return df

def filter_fxsize(fx_df,sc_mean_cutoff=sc_mean_cutoff):
	'''filters for goi'''
	print(fx_df.columns)
	filtered_df=fx_df[fx_df.index.isin(goi)]
	return filtered_df

def makeCompleteName(yName):
	return yName+"\\"+gene_dict[yName]


def transform_and_norm_data(fx_df,fxsize):
	'''drops unneeded columns, normalizes to gen'''
	#drop all but fx size and gene name
	fx_df['one_sided_fx']=fx_df['sp_mean']-fx_df['sc_mean']
	fx_df=fx_df.reset_index()
	fx_df['gene']=fx_df.apply(lambda x: makeCompleteName(x['gene']),axis=1)
	fx_df = fx_df.set_index('gene')
	all_col=fx_df.columns.tolist()
		
	col_to_keep=['one_sided_fx','gene']
	col_to_drop=[col for col in all_col if col not in col_to_keep]
	dropped=fx_df.drop(columns=col_to_drop)
	dropped=dropped.rename(columns={'one_sided_fx':fxsize})
	dropped[fxsize]=dropped[fxsize].div(avg_gen[fxsize])
	#dropped['filler']=1
	#print(dropped)

	#pivot
	#pivoted=dropped.pivot(index='None',columns=['effect_size'],values='effect_size')

	#return pivoted
	return dropped

def plotHeatmap(fx_df):
	print(fx_df)
	fig, ax = plt.subplots(figsize=(8,20))
	mask = fx_df.isnull()
	ax=sns.heatmap(fx_df,cmap=cmap,square=False,mask=mask,robust=True,center=0.0)
	plt.tight_layout()
	plt.tick_params(axis='both', which='major', labelsize=10, labelbottom = False, bottom=False, top = False, labeltop=True)
	plt.show()
	#plt.savefig(figName)

	
def getHeatMapInput(fxsize):
	condition=fxsizes[fxsize]
	fxsize_df=parse_fxsize(fxsize)
	filtered_data=filter_fxsize(fxsize_df,sc_mean_cutoff=sc_mean_cutoff)
	transformed_data=transform_and_norm_data(filtered_data,condition)
	#plotHeatmap(transformed_data)
	return transformed_data

###RUN###

if __name__ == "__main__":

	#build gene dict
	ann_dict=ParseFromGFF(gff)

	gene_dict = {}
	for yName in goi:
		geneName=ann_dict[yName]
		gene_dict[yName]=geneName
	#deal with missing common names
	gene_dict['YIL152W']='VPR1'
	gene_dict['YGL082W']='MIY1'
	gene_dict['YLR422W']='DCK1'
	gene_dict['YJR107W']='LIH1'

	#make heatmap
	fxsize_HeatMapInput=[]
	for fxsize in fxsizes:
		fxsize_HeatMapInput.append(getHeatMapInput(fxsize))

	mergedHeatMapInput=pd.concat(fxsize_HeatMapInput, axis=1, join="outer")
	mergedHeatMapInput.sort_values(by=[fxsizes[sortby]],ascending=False,inplace=True)
	#print(mergedHeatMapInput)
	plotHeatmap(mergedHeatMapInput)

			


