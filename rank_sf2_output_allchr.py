import sys
import numpy as np
import matplotlib.pyplot as plt

'''
rank sweepfinder 2 outputs by gene
assuming file name in form: Bergstrom2014_Spar_1.sf2  
'''
#Parameters
gff='/usr2/people/mabrams/Amended_Genomes/Z1/Z1.gff'
save_prefix='Bergstrom2014_ranked_genes'


filenames = sys.argv[1:]


def ParseFromGFF(gfffile):
    '''
    Parses SGD features flat file
    Input: SGD_features.tab file
    Output: dict of {chrom:{gene:[start,stop]}}
    '''

    ann_dict={}
    
    f = open(gfffile)
    lines=[]
    for line in f:
        row_data=line.split("\t")
        chrom=row_data[0]
        start=int(row_data[3])
        stop=int(row_data[4])
        info=row_data[8].split(";")
        yName=info[0].split('=')[1]
        if yName[0]=='Y':
            if chrom in ann_dict:
                ann_dict[chrom][yName]=[start,stop]
            else:
                ann_dict[chrom]={yName:[start,stop]}
    f.close()
    
    return ann_dict

def ParseSF2(sf2_file):
    '''
    Input: outfile from sweepfinder 2
    Output: 
    '''
    pos_list=[]
    LR_list=[]
        
    #add position and coverage for each base to a dictionary for that chromosome
    f = open(sf2_file)
    next(f)
    for line in f:
        row_data = line.strip().split("\t")
        pos_list.append(float(row_data[0]))
        LR_list.append(float(row_data[1]))

    f.close()

    return pos_list, LR_list

def rankGenes(ann_dict, chrom_dict):
    outfileName=save_prefix+'.txt'
    sf2_dict={}
    gene_LRs=[]
    for chrom in ann_dict:
        pos_list=chrom_dict[chrom][0]
        LR_list=chrom_dict[chrom][1]
        for gene in ann_dict[chrom]:
            LRs=[]
            start=ann_dict[chrom][gene][0]
            stop=ann_dict[chrom][gene][1]
            for i in range(len(pos_list)):
                if start<pos_list[i]<stop:
                    LRs.append(LR_list[i])
                    gene_LR=np.mean(LRs)
                    gene_LRs.append(gene_LR)
            sf2_dict[gene]=[gene_LR, chrom, start, stop]

    #print(sf2_dict['YAL054C'])

    with open(outfileName, 'w') as wf:
        wf.writelines('gene\taverageSF2\n')
        for k in sorted(sf2_dict, key=sf2_dict.get, reverse=True):
            #print(gene)
            gene=k
            avgSF2=str(sf2_dict[gene][0])
            chrom=str(sf2_dict[gene][1])
            start=str(sf2_dict[gene][2])
            stop=str(sf2_dict[gene][3]) 
            wf.writelines(gene+'\t'+avgSF2+'\t'+chrom+'\t'+start+'\t'+stop+'\n')

    return sf2_dict,gene_LRs

def plotHistLRs(gene_LRs):

    fig=plt.figure(figsize=(20,15))
    pl = plt.subplot(111)

    a=np.array(gene_LRs)
    pl.hist(a)
    
    pl.set_xlabel('SF2 LR',fontsize=50)
    pl.set_ylabel('number of genes',fontsize=50)

    plt.title(prefix+' LR plot')
    plt.savefig(save_prefix+'_hist.pdf', bbox_inches='tight', format='pdf', dpi=1000)

    plt.clf()
    return None

##RUN###

ann_dict=ParseFromGFF(gff)
chrom_dict={}

#build chrom dict
for sf2file in filenames:
    #get pos and LR
    pos_list, LR_list=ParseSF2(sf2file)
    #get chromosome names
    prefix=sf2file.split('.')[0]
    chrom=prefix.split('_')[2]
    if len(chrom)==1:
        chrom='chr0'+chrom
    else:
        chrom='chr'+chrom
    chrom_dict[chrom]=[pos_list, LR_list]

sf2_dict, gene_LRs= rankGenes(ann_dict, chrom_dict)
plotHistLRs(gene_LRs)
