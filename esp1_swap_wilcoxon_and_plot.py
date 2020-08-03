from scipy import stats
import numpy as np
#import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd

#PARAMETERS:

plot=True #set to True to plot the data
data_file='ESP1_swaps_raw_data.csv'
significance_cutoff=0.05

strains_to_omit=['407', '413']

def parseData(data_file):
    '''parses csv with the data in it
    returns {date: {temp: {strain: [efficiency1, efficiency2]}
    and a dict of {strain: (background,swap)}'''
    data_dict={}
    strain_dict={}
    f=open(data_file)
    next(f)
    for line in f:
        row_data=line.strip().split(',')
        experimenter=row_data[0]
        date=row_data[1]
        strain=row_data[2]
        background=row_data[3]
        swap=row_data[4]
        od_0=float(row_data[5])
        od_1=float(row_data[6])
        temp=row_data[7]

        eff=od_1 - od_0

        if strain not in strains_to_omit:
            #build data dict
            if date in data_dict:
                if temp in data_dict[date]:
                    if strain in data_dict[date][temp]:
                        data_dict[date][temp][strain].append(eff)
                    else:
                        data_dict[date][temp][strain]=[eff]
                else:
                    data_dict[date][temp]={strain:[eff]}
            else:
                data_dict[date]={temp:{strain:[eff]}}

            #build strain key dict
            if strain not in strain_dict:
                strain_dict[strain]=background,swap

                

    return data_dict, strain_dict

def normalizeData(data_dict):
    '''input: dictionary in form {date: {temp: {strain: [efficiency1, efficiency2]}
        output: dictionary {temp: {strain: norm_efficiency 1, norm_efficieny 2}}
        where normalization is wrt the average sc value for that day'''

    norm_dict={'28C':{},'39C':{}}
    for date in data_dict:
        for temp in data_dict[date]:
            average_sc=np.mean(data_dict[date][temp]['68'])
            for strain in data_dict[date][temp]:
                normalized_eff=[]
                for eff in data_dict[date][temp][strain]:
                    normalized_eff.append(float(eff)/float(average_sc))
                if strain in norm_dict[temp]:
                    norm_dict[temp][strain]+=normalized_eff
                else:
                    norm_dict[temp][strain]=normalized_eff

    #print(norm_dict)
    
    return norm_dict

    
def calcAllWilcoxons(norm_dict,strain_dict,data_file):
    outfile_name=data_file.split('.')[0]+'_wilcoxonResults.txt'
    p_dict={'28C':{},'39C':{}}
    with open(outfile_name,'w') as wf:
        for temp in norm_dict:
            print('==='+temp+'===')
            wf.writelines('==='+temp+'===')
            print('background\tswap\tstrain\tp-value')
            wf.writelines('background\tswap\tstrain\tp-value')
            
            ctr_vals=norm_dict[temp]['68']
            avg_ctr_val=np.mean(ctr_vals)
            for strain in norm_dict[temp]:
                if strain!='68':
                    test_vals=norm_dict[temp][strain]
                    #w,p=stats.wilcoxon(ctr_vals,test_vals)  #trying to do this led to ERROR, UNEQUAL # OF VALS -> do JUST test minus average 68
                    test_diffs_from_avg_68=[]
                    for measurement in test_vals:
                        test_diffs_from_avg_68.append(measurement - avg_ctr_val)
                    w,p=stats.wilcoxon(test_diffs_from_avg_68, alternative='less')
                    background=strain_dict[strain][0]
                    swap=strain_dict[strain][1]
                    print(background+'\t'+swap+'\t'+strain+'\t'+ str(p))
                    wf.writelines(background+'\t'+swap+'\t'+strain+'\t'+ str(p))
                    p_dict[temp][strain]=p


##                    if temp=='28C' and strain=='62': #test 62
##                        print(test_vals)
##                        print(test_diffs_from_avg_68)

    return(p_dict)

def plotOnePlot(temp, temp_data,significant_strains):

    w=0.8 #bar width
    colors = [(1, 0.5, 0, 1),(0, 0.5, 1, 1), (0.05, 0.8, 0.1, 1),(0.05, 0.8, 0.1, 1), (0.05, 0.8, 0.1, 1),(0.05, 0.8, 0.1, 1), (0.05, 0.8, 0.1, 1)]
    strains=list(temp_data.keys())
    x=range(len(strains)) #x: num strains
    y=[] #y: mean of bioreps for the strains
    labels= [ ] #label: background, swap
    for strain in strains:
        #y.append(np.mean(temp_data[strain]))
        y.append(np.array(temp_data[strain]))
        labels.append(strain_dict[strain][0]+','+strain_dict[strain][1])
    #print(y)
    fig, ax = plt.subplots()
    ax.bar(x,
       height=[np.mean(yi) for yi in y],
       width=w,    # bar width
       tick_label=labels, color=(0,0,0,0),  # face color transparent
       edgecolor=colors,
       )


    for i in range(len(x)):
        # distribute scatter randomly across whole width of bar
        ax.scatter(x[i] + np.random.random(y[i].size) * w - w / 2, y[i], color=colors[i])


    x_for_line_at_1 = [ax.patches[0].get_x(), ax.patches[-1].get_x() + ax.patches[-1].get_width()] #https://stackoverflow.com/questions/38017465/how-to-add-a-line-on-top-of-a-bar-chart
    ax.plot(x_for_line_at_1, [1,1], 'r--', c='grey', linewidth=1)

    plt.title(temp,fontsize=20,loc='center')
    plt.ylabel('Normalized Efficiency', fontsize=20)
    plt.xlabel('Strain (Background, Swap)', fontsize=20)

    plt.savefig(temp+"_esp1_swaps_barplot.png")
    
##    np.random.seed(123) #bars with scatter ex from https://stackoverflow.com/questions/51027717/pyplot-bar-charts-with-individual-data-points/51032760
##
##    w = 0.8    # bar width
##    x = [1, 2] # x-coordinates of your bars
##    colors = [(0, 0, 1, 1), (1, 0, 0, 1)]    # corresponding colors
##    y = [np.random.random(30) * 2 + 5,       # data series
##        np.random.random(10) * 3 + 8]
##
##    print(y)
##
##    fig, ax = plt.subplots()
##    ax.bar(x,
##           height=[np.mean(yi) for yi in y],
##           yerr=[np.std(yi) for yi in y],    # error bars
##           capsize=12, # error bar cap width in points
##           width=w,    # bar width
##           tick_label=["control", "test"],
##           color=(0,0,0,0),  # face color transparent
##           edgecolor=colors,
##           #ecolor=colors,    # error bar colors; setting this raises an error for whatever reason.
##           )
##
##    for i in range(len(x)):
##        # distribute scatter randomly across whole width of bar
##        ax.scatter(x[i] + np.random.random(y[i].size) * w - w / 2, y[i], color=colors[i])

        

    

    
    return                
                
def plotBars(norm_dict,strain_dict,data_file,p_dict,significance=0.05):
    for temp in norm_dict:
        temp_data=norm_dict[temp]
        significant_strains=[]
        for strain in p_dict[temp]:
            if p_dict[temp][strain]<significance:
                significant_strains.append(strain)
        plotOnePlot(temp, temp_data, significant_strains)
    return
        
        

        
###START###
print('done loading libraries')
data_dict,strain_dict=parseData(data_file)
print('done parsing data')
norm_dict=normalizeData(data_dict)
print('done normalizing data')
p_dict=calcAllWilcoxons(norm_dict,strain_dict,data_file)
if plot==True:
    plotBars(norm_dict, strain_dict, data_file,p_dict,significance=significance_cutoff)
    print('done plotting')
