import sys
import numpy as np


#usage: python3 10effect_size.py rep_data.insert_ratios

#purpose: calculate effect size (distance between means) of sp and sc observations

#Parameters

p_cutoff=0.05 #p-value cutoff to analyze effect size.
sp_max_offcenter=0.75 #cutoff for how far the sp distribution can be from zero (mean across ALL observations).  
save_fx_size=True #save the effect sizes of all genes analyzed in an outfile




def ParseFitness(rep_data):
    '''get the replicates' fitness data'''
    sp_dict={}
    sc_dict={}


    print('...parsing logratios...')

    f = open(rep_data)
    header_line=True
    for line in f:
        row_data=line.split("\t")
        
        if header_line== True:
            annotation_index=row_data.index('annotation')
            logr_indices = []
            for item in row_data:
                if '39_28_log2' in item and 'br' not in item:  #take individual bioreps' logratios, not the logratio of averaged brs or its cv
                    logr_indices.append(row_data.index(item))
            header_line = False
            
        else:
            annotation=row_data[annotation_index]
            allele=annotation[:2]
            gene=annotation[2:]
            for logr_index in logr_indices:
                logr=row_data[logr_index]
                #add to sc or sp dict
                if logr!='':
                    if allele=='sp':
                        if gene in sp_dict:
                            sp_dict[gene].append(float(logr))
                        else:
                            sp_dict[gene]=[float(logr)]
                    elif allele=='sc':
                        if gene in sc_dict:
                            sc_dict[gene].append(float(logr))
                        else:
                            sc_dict[gene]=[float(logr)]
            
                   
    f.close()

    print('...done parsing logratios...')

##    print('YLR397C sp_ratios')
##    print(sp_dict['YLR397C'])


    return sp_dict,sc_dict


def CalcEffectSize(sc_dict,sp_dict,save_fx_size=False,save_name='effect_sizes.txt',printTop=5):
    '''calculates effect size for all significant genes'''

    print('...calculating effect sizes...')
    fx_dict={}
    mean_dict={}
    for gene in sc_dict:
        if gene in sp_dict:
            sc_mean=np.mean(sc_dict[gene])
            sp_mean=np.mean(sp_dict[gene])
            if gene=='YLR397C':
                        print('sp_mean')
                        print(sp_mean)
            fx_dict[gene]=abs(sc_mean-sp_mean) #fx size is absolute value of difference between means of individual brs logratios
            mean_dict[gene]=[sc_mean,sp_mean]

    print('...done calculating effect sizes!')
    
    if printTop>0: #print the top some number of effect sizes
        print('top '+str(printTop)+' effect sizes:')
        for i in range(printTop):
            for gene in sorted(fx_dict, key=fx_dict.get, reverse=True):
              if printTop>0:
                  print(gene, fx_dict[gene])
                  printTop-=1
    if save_fx_size==True:
        print('...saving outfile:')
        wf=open(save_name,'w')
        wf.writelines('gene\teffect_size\tsc_mean\tsp_mean\tmwu_pval\n')
        for gene in sorted(fx_dict, key=fx_dict.get, reverse=True):
            wf.writelines(gene +'\t'+str(fx_dict[gene])+'\t'+str(mean_dict[gene][0])+'\t'+str(mean_dict[gene][1])+'\n')
        wf.close()
    return fx_dict,mean_dict



###RUN###

if __name__ == "__main__":

    rep_data=sys.argv[1]
 
    sp_dict,sc_dict=ParseFitness(rep_data)

    effect_dict,mean_dict=CalcEffectSize(sc_dict,sp_dict,save_fx_size=True,save_name=rep_data.split('.insert_ratios')[0]+'.all_fxsizes',printTop=0)
 

