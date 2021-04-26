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

##USAGE## python heatmap.py  #make a simple heatmap (one effect size per condition)


##PARAMETERS##

#fxsize files:
fx36='36C_2.0_1.1_coeffvarANDreadcutoff_uncollapsedMBA_unpaired.all_fxsizes' #all 36C data
fx37='37C_2.0_1.1_coeffvarANDreadcutoff_uncollapsedMBA_unpaired.all_fxsizes' #all 37C data 


sortby=fx37
#fxsizes={fx36:'36C',fx37:'37C',fx38:'38C'}
fxsizes={fx36:'36C',fx37:'37C'}

gen_36=[12.84584039,12.96574443,12.82881066,12.47694361,12.69612127,12.57667477,12.72712914,12.42038924,12.4493805,12.52735329,13.08141686,13.31136996]
gen_37=[10.71381697,10.6059797,10.70591248,10.74969728,10.72307974,10.87017294,10.95210186,10.61683005,9.636340693,10.18646952,10.72492607,10.55877579]

avg_gen={'36C':np.mean(gen_36),'37C':np.mean(gen_37)} 

goi=['YLR397C','YGR098C','YKR054C','YHR023W',
     'YDR180W','YCR042C','YNL172W'] #thermotolerance loci from Weiss et al 2018, excluding YMR168C b/c only one br


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
	return df

def filter_fxsize(fx_df):
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
	return dropped

def plotHeatmap(fx_df):
	print(fx_df)
	fig, ax = plt.subplots(figsize=(8,20))
	mask = fx_df.isnull()
	ax=sns.heatmap(fx_df,cmap=cmap,square=False,mask=mask,robust=True,center=0.0)
	plt.tick_params(axis='both', which='major', labelsize=10, labelbottom = False, bottom=False, top = False, labeltop=True,rotation=0)
	#plt.tight_layout()
	plt.show()
	#plt.savefig(figName)

	
def getHeatMapInput(fxsize):
	condition=fxsizes[fxsize]
	fxsize_df=parse_fxsize(fxsize)
	filtered_data=filter_fxsize(fxsize_df)
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
		
	fxsize_HeatMapInput=[]
	for fxsize in fxsizes:
		fxsize_HeatMapInput.append(getHeatMapInput(fxsize))

	mergedHeatMapInput=pd.concat(fxsize_HeatMapInput, axis=1, join="outer")
	mergedHeatMapInput.sort_values(by=[fxsizes[sortby]],ascending=False,inplace=True)
	#print(mergedHeatMapInput)
	plotHeatmap(mergedHeatMapInput)

			


