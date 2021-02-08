import sys
import numpy as np


#usage: python3 10effect_size.py rep_data.txt y.test_test_results

#purpose: calculate effect size (distance between means) of sp and sc observations

#Parameters

p_cutoff=0.05 #p-value cutoff to analyze effect size.
sp_max_offcenter=0.75 #cutoff for how far the sp distribution can be from zero
save_fx_size=True #save the effect sizes of all genes analyzed in an outfile

def ParseStats(stats_data, cutoff=0.05):
    '''get dict of pvalues of genes below significance cutoff'''

    print('...parsing stats...')
    stats_dict={}
    with open (stats_data) as f:
        next(f) #skip header line
        for line in f:
            row_data=line.split("\t")
            gene = row_data[0]
            adj_pval=float(row_data[2].strip('\n'))
            stats_dict[gene]=adj_pval

    print('...done parsing stats...')
    return stats_dict


def ParseFitness(rep_data,stats_dict):
    '''get the replicates' fitness data'''
    sp_dict={}
    sc_dict={}

    for key in stats_dict:
        sp_dict[key]=[]
        sc_dict[key]=[]

    print('...parsing logratios...')

    f = open(rep_data)
    header_line=True
    for line in f:
        row_data=line.split("\t")
        
        if header_line== True:
            gene_index = row_data.index('gene')
            allele_index = row_data.index('allele')
            logr_indices = []
            for item in row_data:
                if '39_28_log2' in item and 'br_averaged_reads' not in item: #take individual bioreps' logratios, not the logratio of averaged brs
                    logr_indices.append(row_data.index(item))
            header_line = False
            
        else:
            gene=row_data[gene_index]
            if gene in stats_dict:
                allele=row_data[allele_index]
                for logr_index in logr_indices:
                    logr=row_data[logr_index]
                    #add to sc or sp dict
                    if logr!='':
                        if allele=='sp':
                            sp_dict[gene].append(float(logr))
                        elif allele=='sc':
                            sc_dict[gene].append(float(logr))
    f.close()

    print('...done parsing logratios...')

    return sp_dict,sc_dict


def CalcEffectSize(stats_dict,sc_dict,sp_dict,save_fx_size=False,save_name='effect_sizes.txt',printTop=5):
    '''calculates effect size for all significant genes'''

    print('...calculating effect sizes...')
    fx_dict={}
    for gene in stats_dict:
        sc_mean=np.mean(sc_dict[gene])
        sp_mean=np.mean(sp_dict[gene])
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
            wf.writelines(gene +'\t'+str(fx_dict[gene])+'\t'+str(mean_dict[gene][0])+'\t'+str(mean_dict[gene][1])+'\t'+str(stats_dict[gene])'\n')
        wf.close()
    return fx_dict,mean_dict


def FilterEffectSize(effect_dict, mean_dict,stats_dict,p_cutoff,sp_max_min,save_name='filtered_effect_sizes.txt',printTop=4):
    filtered_effect_dict={}
    
    significant=0 #filter by p
    outside_sp_offcenter=0
    for gene in stats_dict:
        if stats_dict[gene]>p_cutoff: #filter by p
            significant+=1
            if abs(mean_dict[gene][2]<sp_min_max: #filter by how off-center sp is
                   filtered_effect_dict[gene]=effect_dict[gene]
            else:
                outside_sp_offcenter+=1

    print(str(significant)+' genes significant at p-threshold of '+str(p_cutoff))
    print(str(outside_sp_offcenter)+' hits discarded out because sp too offcenter (abs sp mean >'+str(sp_min_max)+')')

    sc_greater=0
    sp_greater=0
    for gene in filtered_effect_dict:
        if mean_dict[gene][0]>mean_dict[gene][1]:
            sc_greater+=1
        else:
            sp_greater+=1
            
    print(str(sc_greater)+' filtered genes where disruptions in sc had a higher fitness')
    print(str(sp_greater)+' filtered genes where disruptions in sp had a higher fitness')

                
    if printTop>0: #print the top some number of effect sizes
        print('top '+str(printTop)+'filtered effect sizes:')
        for i in range(printTop):
            for gene in sorted(filtered_effect_dict, key=filtered_effect_dict.get, reverse=True):
              if printTop>0:
                  print(gene, filtered_effect_dict[gene])
                  printTop-=1
    if save_fx_size==True:
        print('...saving outfile:')
        wf=open(save_name,'w')
        wf.writelines('gene\teffect_size\tsc_mean\tsp_mean\tmwu_pval\n')
        for gene in sorted(filtered_effect_dict, key=filtered_effect_dict.get, reverse=True):
            wf.writelines(gene +'\t'+str(filtered_effect_dict[gene])+'\t'+str(mean_dict[gene][0])+'\t'+str(mean_dict[gene][1])+'\t'+str(stats_dict[gene])'\n')
        wf.close()

    
    return filtered_effect_dict
        

##def PlotDensities(sp_dict,sc_dict,stats_dict=None,specific_genes=None):
##    if stats_dict!=None:
##        for gene in stats_dict:
##            PlotDensity(sp_dict[gene], sc_dict[gene],gene)
##    elif specific_genes!=None:
##        for gene in specific_genes:
##            if gene in sp_dict and gene in sc_dict:
##                PlotDensity(sp_dict[gene], sc_dict[gene],gene)
##            else:
##                print(gene+" is not in both alleles")
##    return 'done plotting'

###RUN###

if __name__ == "__main__":

    rep_data=sys.argv[1]
    stats_data=sys.argv[2]

    stats_dict=ParseStats(stats_data)
    sp_dict,sc_dict=ParseFitness(rep_data,stats_dict)

    effect_dict,mean_dict=CalcEffectSize(stats_dict,sc_dict,sp_dict,save_fx_size=True,save_name=stats_data.split('.mwu')[0]+'.all_fxsizes',printTop=0)
    filtered_effect_dict=FilterEffectSize(effect_dict, mean_dict,stats_dict,p_cutoff,sp_max_min,save_name=stats_data.split('.mwu')[0]+'p'+str(p_cutoff)+'spmaxoffby'+str(sp_min_max).fxsizes',printTop=4)
    #PlotDensities(sp_dict,sc_dict,stats_dict=stats_dict,specific_genes=specific_genes_to_include)

