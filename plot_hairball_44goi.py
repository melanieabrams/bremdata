import networkx as nx
import requests
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
# from matplotlib import cm
from nxviz.plots import CircosPlot


import matplotlib.patches as mpatches
from matplotlib.colors import to_hex
csfont = {'fontname':'Arial'}
plt.rcParams.update({
    'font.sans-serif': 'Arial',
    'font.family': 'sans-serif',
    'font.size':10
})
#tutorial: https://towardsdatascience.com/chord-diagrams-of-protein-interaction-networks-in-python-9589affc8b91
#tutorial written by Ford Combs
#adding legend from https://github.com/ericmjl/nxviz/issues/592
#legend code from Walter Hernandez

##Parameters##
goi=['YGR198W', 'YMR207C', 'YGL082W', 'YNL049C', 'YDL035C', 'YDR508C', 'YBR136W',
 'YML099C', 'YPL254W', 'YIL152W', 'YKL017C', 'YGR140W', 'YJR127C', 'YDR375C',
 'YOR091W', 'YLR397C', 'YNL132W', 'YMR078C', 'YLR422W', 'YMR125W', 'YOR371C',
 'YMR094W', 'YMR167W', 'YDR103W', 'YDR318W', 'YAL026C', 'YDR180W', 'YOR092W',
 'YDR235W', 'YER151C', 'YMR275C', 'YKL114C', 'YOL081W', 'YPR049C', 'YGL095C',
 'YDR456W', 'YKL197C', 'YIL068C', 'YOR326W', 'YNR045W', 'YJR107W', 'YPL268W',
 'YJL062W', 'YCR042C'] #hits 37_2.0_1.1_..._10.0_3.0 with sc defect and effect size >0.5

gff='saccharomyces_cerevisiae_R64-1-1_20110208_withoutFasta.gff'


##functions##

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

def makeCompleteName(geneName):
	return gene_dict[geneName]+"\\"+'\n'+geneName

##Run##


#build gene dict
ann_dict=ParseFromGFF(gff)

gene_dict = {}
for yName in goi:
        geneName=ann_dict[yName]
        gene_dict[geneName]=yName
#deal with missing common names
gene_dict['VPR1']='YIL152W'
gene_dict['MIY1']='YGL082W'
gene_dict['DCK1']='YLR422W'
gene_dict['LIH1']='YJR107W'


#Gather the data

#use python request to get interactions from STRING db

proteins = '%0d'.join(goi)
url = 'https://string-db.org/api/tsv/network?identifiers=' + proteins + '&species=4932'
#print(url)

r = requests.get(url)


## build a data frame in which each row contains a pair of proteins and the interaction score between them.

lines = r.text.split('\n') # pull the text from the response object and split based on new lines
data = [l.split('\t') for l in lines] # split each line into its components based on tabs
# convert to dataframe using the first row as the column names; drop empty, final row
df = pd.DataFrame(data[1:-1], columns = data[0]) 
# dataframe with the preferred names of the two proteins and the score of the interaction
interactions = df[['preferredName_A', 'preferredName_B', 'score']]
interactions['preferredName_A']=interactions.apply(lambda x: makeCompleteName(x['preferredName_A']),axis=1)
interactions['preferredName_B']=interactions.apply(lambda x: makeCompleteName(x['preferredName_B']),axis=1)
print(interactions)

#build the graph using NetworkX

G=nx.Graph(name='Protein Interaction Graph')
interactions = np.array(interactions) # convert to array for clarity
for i in range(len(interactions)):
    interaction = interactions[i]
    a = interaction[0] # protein a node
    b = interaction[1] # protein b node
    w = int(float(interaction[2])*100) # score as weighted edge
    
    # To include all the weighted connections, uncomment the following line
    G.add_weighted_edges_from([(a,b,w)])
    
##    # To only keep high scoring edges, use the following lines
##    if w > 80: # only keep high scoring edges
##        G.add_weighted_edges_from([(a,b,w)])


###Draw a simple  plot:
##c = CircosPlot(G,node_labels=True)
##c.draw()
##plt.show()


#Improved plot

# function to rescale list of values to range [newmin,newmax]
def rescale(l,newmin,newmax,rnd=False):
    arr = list(l)
    return [round((x-min(arr))/(max(arr)-min(arr))*(newmax-newmin)+newmin,2) for x in arr]

nodelist = [n for n in G.nodes]
ws = rescale([float(G[u][v]['weight']) for u,v in G.edges],1,10)
# alternative method below (uncommented b/c kept all weight)
ws = rescale([float(G[u][v]['weight'])**70 for u,v in G.edges],1,50)
edgelist = [(str(u),str(v),{"weight":ws.pop(0)}) for u,v in G.edges]

# create new graph using nodelist and edgelist
g = nx.Graph(name='Protein Interaction Graph')
g.add_nodes_from(nodelist)
g.add_edges_from(edgelist)
# go through nodes in graph G and store their degree as "class" in graph g
for v in G:
    g.nodes[v]["class"] = G.degree(v)

#draw
c = CircosPlot(graph=g,figsize=(20,20),node_grouping="class", node_color="class",
               edge_width="weight",node_labels=True,fontsize=10,fontfamily='sans-serif')#,
##               group_label_offset=0.75,group_label_position='beginning',group_legend=True,
##               group_label_color=True)
c.figure.tight_layout() 
c.draw()


 #Get the label and color for each group used by nxviz
seen = set()
colors_group = [x for x in c.node_colors if not (x in seen or seen.add(x))] #Gets colors in RGBA
labels_group = sorted(list(set([g.nodes[n][c.node_color] for n in G.nodes])))

#Create patchList to use as handle for plt.legend()
patchList = []
for color, label in zip(colors_group, labels_group):
    color = to_hex(color, keep_alpha=True) #Convert RGBA to HEX value
    data_key = mpatches.Patch(color=color, label=label)
    patchList.append(data_key)

#Set the labels with the custom patchList
plt.legend(handles=patchList,
           loc="lower center",
           title="Interaction Class",
           ncol=4, # 4 columns to spread the legends
           borderpad=1,
           bbox_to_anchor =(1.05, -0.0), # To move the legend box lower than the graph
           shadow=False,
           fancybox=True)

plt.show()
