import sys
import numpy as np
from funcy import flatten, isa

#usage: python3 10effect_size.py rep_data.insert_ratios

#purpose: calculate effect size (distance between means) of sp and sc observations

##DECIDED THIS IS NOT GREAT BECAUSE IT EITHER COLLAPSES BY INS, BR, OR HAS A TON OF NAN.

#Parameters

p_cutoff=0.05 #p-value cutoff to analyze effect size.
sp_max_offcenter=0.75 #cutoff for how far the sp distribution can be from zero (mean across ALL observations for sp, compared to ONE observation for sc).  
save_fx_size=True #save the effect sizes of all genes analyzed in an outfile
temp='37C'




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
            if gene not in sp_dict:
                sp_dict[gene]=[]
                for i in range(len(logr_indices)):
                    sp_dict[gene].append([])
            if gene not in sc_dict:
                sc_dict[gene]=[]
                for i in range(len(logr_indices)):
                    sc_dict[gene].append([])
            for i in range(len(logr_indices)):
                logr_index=logr_indices[i]
                logr=row_data[logr_index]
                #add to sc or sp dict           
                if logr!='':
                    if allele=='sp':
                        if gene in sp_dict:
                            sp_dict[gene][i].append(float(logr))
                        else:
                            sp_dict[gene][i]=[float(logr)]
                    elif allele=='sc':
                        if gene in sc_dict:
                            sc_dict[gene][i].append(float(logr))
                        else:
                            sc_dict[gene][i]=[float(logr)]
                else:
                    if allele=='sp':
                        if gene in sp_dict:
                            sp_dict[gene][i].append(np.nan)
                        else:
                            sp_dict[gene][i]=[np.nan]
                    elif allele=='sc':
                        if gene in sc_dict:
                            sc_dict[gene][i].append(np.nan)
                        else:
                            sc_dict[gene][i]=[np.nan]
            
                   
    f.close()

    print('...done parsing logratios...')

##    print('YLR397C sp_ratios')
##    print(sp_dict['YLR397C'])


    return sp_dict,sc_dict,logr_indices


def CalcEffectSize(sc_dict,sp_dict,logr_indices,save_fx_size=False,save_name='effect_sizes.txt',printTop=5):
    '''calculates effect size for all significant genes'''

    print('...calculating effect sizes...')
    fx_dict={}
    mean_dict={}
    for gene in sc_dict:
        if gene in sp_dict:
            if len(sc_dict[gene])>0 and len(sp_dict[gene])>0:
                all_sp_obs=list(flatten(sp_dict[gene]))
                all_sc_obs=list(flatten(sc_dict[gene]))
                sp_mean=np.nanmean(all_sp_obs)
                sc_mean=np.nanmean(all_sc_obs)
                fx_dict[gene]=[]
                
                for i in range(len(logr_indices)): #fx size is absolute value of difference between means of individual brs logratios
                    if sc_dict[gene][i]==np.nan:
                        fx_dict[gene].append(np.nan)
                    else:
                        fx_dict[gene].append(sp_mean-np.nanmean(sc_dict[gene][i]))
                fx_dict[gene].append(sp_mean-sc_mean) #also add mean for context https://github.com/mwaskom/seaborn/issues/375

    print('...done calculating effect sizes!')

    print(sc_dict['YLR397C'])
    
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
        header='gene'
        for i in range(len(logr_indices)):
            header+='\t'+temp+'_'+str(i+1)
        wf.writelines(header+'\tmean\n')
        
        for gene in sorted(fx_dict, key=fx_dict.get, reverse=True):
            gene_str=gene
            for i in range(len(fx_dict[gene])):
                gene_str+='\t'+str(fx_dict[gene][i])
            wf.writelines(gene_str+'\n')
        wf.close()
    return fx_dict,mean_dict



###RUN###

if __name__ == "__main__":

    rep_data=sys.argv[1]
 
    sp_dict,sc_dict,logr_indices=ParseFitness(rep_data)

    effect_dict,mean_dict=CalcEffectSize(sc_dict,sp_dict,logr_indices,save_fx_size=True,save_name=temp+'_'+rep_data.split('.insert_ratios')[0]+'.br_fxsizes',printTop=0)
 

