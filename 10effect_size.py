import sys
import numpy as np


#usage: python3 10effect_size.py rep_data.txt y.test_test_results

#purpose: calculate effect size (distance between means) of sp and sc observations

#Parameters

cutoff=0.05 #p-value cutoff to analyze effect size.  pass '1' here to analyze all genes' effect sizes.
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
            if adj_pval<cutoff:
                stats_dict[gene]=adj_pval

    print('...done parsing stats...')
    return stats_dict


def ParseFitness(rep_data,stats_dict):
    '''get the replicates' fitness data for all genes below significance cutoff'''
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


def CalcEffectSize(stats_data,sc_dict,sp_dict,save_fx_size=False,save_name='effect_sizes.txt',printTop=5):
    '''calculates effect size for all significant genes'''

    print('...calculating effect sizes...')
    fx_dict={}
    for gene in stats_dict:
        sc_mean=np.mean(sc_dict[gene])
        sp_mean=np.mean(sp_dict[gene])
        fx_dict[gene]=abs(sc_mean-sp_mean) #fx size is absolute value of difference between means of individual brs logratios

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
        wf.writelines('gene\teffect_size\n')
        for gene in sorted(fx_dict, key=fx_dict.get, reverse=True):
            wf.writelines(gene +'\t'+str(fx_dict[gene])+'\n')
        wf.close()
    return fx_dict
        
        

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
    
    stats_dict=ParseStats(stats_data,cutoff=cutoff)
    sp_dict,sc_dict=ParseFitness(rep_data,stats_dict)

    effect_dict=CalcEffectSize(stats_data,sc_dict,sp_dict,save_fx_size=True,save_name=stats_data.split('.mwu')[0]+'.fxsizes',printTop=4)
    #PlotDensities(sp_dict,sc_dict,stats_dict=stats_dict,specific_genes=specific_genes_to_include)

