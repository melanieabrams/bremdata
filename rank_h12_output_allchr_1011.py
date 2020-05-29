import sys
import numpy as np
import matplotlib.pyplot as plt

'''
rank selecthapstats h12 outputs by gene
assuming file name in form: 
merged_6AfricanBeer_chromosome1.h12

'''
#Parameters
#gff='/usr2/people/mabrams/Amended_Genomes/D1373/DBVPG1373.gff'
gff='/usr2/people/mabrams/Amended_Genomes/S288C/saccharomyces_cerevisiae_R64-1-1_20110208_withoutFasta.gff'
roman_numerals_in_gff=True



filenames = sys.argv[1:]
save_prefix=filenames[0].split('_chromosome')[0]


def ParseFromGFF(gfffile):
    '''
    Parses SGD features flat file
    Input: SGD_features.tab file
    Output: dict of {chrom:{gene:[start,stop]}}
    '''
    roman_to_numerals={
        'chrI':'chr01','chrII':'chr02','chrIII':'chr03','chrIV':'chr04','chrV':'chr05',
        'chrVI':'chr06','chrVII':'chr07','chrVIII':'chr08','chrIX':'chr09','chrX':'chr10',
        'chrXI':'chr11','chrXII':'chr12','chrXIII':'chr13','chrXIV':'chr14','chrXV':'chr15',
        'chrXVI':'chr16','chrMito':'chrMito','2-micron':'chr2u'}
    ann_dict={}
    
    f = open(gfffile)
    lines=[]
    for line in f:
        if line[0]!='#': #skip header rows
            row_data=line.split("\t")
            #print(row_data)
            if roman_numerals_in_gff==True: #convert roman numerals if needed
                chrom=roman_to_numerals[row_data[0]]
            else:
                chrom=row_data[0]
            start=int(row_data[3])
            stop=int(row_data[4])
            info=row_data[8].split(";")
            yName=info[0].split('=')[1]
            #print(yName)
            if yName[0]=='Y' and len(yName)>5:
                if chrom in ann_dict:
                    ann_dict[chrom][yName]=[start,stop]
                else:
                    ann_dict[chrom]={yName:[start,stop]}
    f.close()
    
    return ann_dict

def ParseH12(h12_file):
    '''
    Input: outfile from sweepfinder 2
    Output: 
    '''
    pos_list=[]
    H12_list=[]
    ratioH2H1_list=[]
        
    #add position and values for each base to a dictionary for that chromosome
    f = open(h12_file)
    next(f)
    for line in f:
        row_data = line.strip().split("\t")
        peak_ctr=int(row_data[0])
        pos_list.append(peak_ctr)
        H12_list.append(float(row_data[8]))
        ratioH2H1_list.append(float(row_data[9]))
        #columns:
            #1ctrcoord (index 0)
            #2leftcoord
            #3rightcoord
            #4K
            #55hapfreqspec
            #6strainnum
            #7H1
            #8H2
            #9H12 (index 8)
            #10H2/H1 (index 9)
            #11H123 (index 10)
 
    f.close()

    return pos_list, H12_list, ratioH2H1_list

   

def rankGenes(ann_dict, chrom_dict):
    outfileName=save_prefix+'_ranked.txt'
    H12_dict={}
    gene_H12s=[]
    gene_H2H1ratios=[]
    for chrom in ann_dict:
        pos_list=chrom_dict[chrom][0]
        H12_list=chrom_dict[chrom][1]
        ratio_H2H1_list=chrom_dict[chrom][2]
        for gene in ann_dict[chrom]:
            H12s=[]
            H2H1_ratios=[]
            start=ann_dict[chrom][gene][0]
            stop=ann_dict[chrom][gene][1]
            for i in range(len(pos_list)):
                if start<pos_list[i]<stop:
##                    if gene=="YAL056W":
##                        print(gene)
##                        print("H12 list")
##                        print(H12_list[i])
##                        print("H2/H1 ratios")
##                        print(ratio_H2H1_list[i])
                    H12s.append(H12_list[i])
                    H2H1_ratios.append(ratio_H2H1_list[i])
                    

            if len(H12s)>0:
                gene_H12=np.mean(H12s)
                gene_H12s.append(gene_H12)
                gene_H2H1ratio=np.mean(H2H1_ratios)
                gene_H2H1ratios.append(gene_H2H1ratio)
                H12_dict[gene]=[gene_H12, gene_H2H1ratio, chrom, start, stop]

    #print(sf2_dict['YAL054C'])

    with open(outfileName, 'w') as wf:
        wf.writelines('gene\taverageH12\taverageH2H1ratio\tchrom\tstart\tstop\n')
        for k in sorted(H12_dict, key=H12_dict.get, reverse=True):
            #print(gene)
            gene=k
            avgH12=str(H12_dict[gene][0])
            avgH2H1ratio=str(H12_dict[gene][1])
            chrom=str(H12_dict[gene][2])
            start=str(H12_dict[gene][3])
            stop=str(H12_dict[gene][4]) 
            wf.writelines(gene+'\t'+avgH12+'\t'+avgH2H1ratio+'\t'+chrom+'\t'+start+'\t'+stop+'\n')

    return H12_dict,gene_H12s, gene_H2H1ratios

def plotHistHstats(gene_stats, stat):

    fig=plt.figure(figsize=(20,15))
    pl = plt.subplot(111)

    a=np.array(gene_stats)
    pl.hist(a)
    
    pl.set_xlabel('SelectHapStats '+stat,fontsize=50)
    pl.set_ylabel('number of genes',fontsize=50)

    plt.title(prefix+' '+stat+' histogram')
    plt.savefig(save_prefix+stat+'_hist.pdf', bbox_inches='tight', format='pdf', dpi=1000)

    plt.clf()
    return None



##RUN###

ann_dict=ParseFromGFF(gff)
chrom_dict={}

#build chrom dict
for h12file in filenames:
    #get pos and vals
    pos_list, H12_list, ratioH2H1_list=ParseH12(h12file)
    #get chromosome names
    prefix=h12file.split('.')[0]
    chrom=prefix.split('chromosome')[1]
    if len(chrom)==1:
        chrom='chr0'+chrom
    else:
        chrom='chr'+chrom
    chrom_dict[chrom]=[pos_list, H12_list, ratioH2H1_list]

h12_dict, gene_h12s, gene_H2H1ratios= rankGenes(ann_dict, chrom_dict)
plotHistHstats(gene_h12s, 'H12')
plotHistHstats(gene_H2H1ratios, 'H2H1ratio')
