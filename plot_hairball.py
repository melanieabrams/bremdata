import networkx as nx
import requests
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
# from matplotlib import cm
from nxviz.plots import CircosPlot


import matplotlib.patches as mpatches
from matplotlib.colors import to_hex

#tutorial: https://towardsdatascience.com/chord-diagrams-of-protein-interaction-networks-in-python-9589affc8b91
#tutorial written by Ford Combs
#adding legend from https://github.com/ericmjl/nxviz/issues/592
#legend code from Walter Hernandez

##Parameters##
goi=['YGR198W','YLR397C','YBR136W','YGR140W','YNL049C','YKL017C','YGL082W',
     'YMR125W','YDR508C','YMR207C','YDR375C','YDR180W','YKL197C','YDR318W',
     'YMR094W','YOR326W','YBR081C','YPR049C','YIL152W','YER151C','YJR107W',
     'YAL026C','YDR456W','YLR141W','YPL268W','YDR235W']

##Run##

#Gather the data

#use python request to get interactions from STRING db

proteins = '%0d'.join(goi)
url = 'https://string-db.org/api/tsv/network?identifiers=' + proteins + '&species=4932'
print(url)

r = requests.get(url)


## build a data frame in which each row contains a pair of proteins and the interaction score between them.

lines = r.text.split('\n') # pull the text from the response object and split based on new lines
data = [l.split('\t') for l in lines] # split each line into its components based on tabs
# convert to dataframe using the first row as the column names; drop empty, final row
df = pd.DataFrame(data[1:-1], columns = data[0]) 
# dataframe with the preferred names of the two proteins and the score of the interaction
interactions = df[['preferredName_A', 'preferredName_B', 'score']]

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
c = CircosPlot(graph=g,figsize=(13,13),node_grouping="class", node_color="class",
               edge_width="weight",node_labels=True,fontsize=11,fontfamily='sans-serif')#,
##               group_label_offset=0.75,group_label_position='beginning',group_legend=True,
##               group_label_color=True)
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
