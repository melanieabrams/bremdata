import sys
import numpy as np
import matplotlib.pyplot as plt

'''
rank selecthapstats g1 (or g12) outputs by gene
assuming file name in form: 
merged_6AfricanBeer_chromosome1.h12_h2h1

'''
#Parameters
gff='/usr2/people/mabrams/Amended_Genomes/SaccSensuStricto/Spar.gff'
roman_numerals_in_gff=False



filenames = sys.argv[1:]
save_prefix=filenames[0].split('_Spar')[0]


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
            SGD=info[4]
            if SGD!='SGD=':
                yName=SGD.split('=')[1]

                if yName[0]=='Y' and len(yName)>5:
                    if chrom in ann_dict:
                        ann_dict[chrom][yName]=[start,stop]
                    else:
                        ann_dict[chrom]={yName:[start,stop]}
    f.close()

    print(ann_dict.keys())
    
    return ann_dict

def ParseH12(h12_file):
    '''
    Input: outfile from sweepfinder 2
    Output: 
    '''
    pos_list=[]
    G1_list=[]

    #add position and values for each base to a dictionary for that chromosome
    f = open(h12_file)
    next(f)
    for line in f:
        row_data = line.strip().split("\t")
        peak_ctr=int(row_data[0])
        pos_list.append(peak_ctr)
        G1_list.append(float(row_data[6]))
        
        #columns:
            #1ctrcoord (index 0)
            #2leftcoord
            #3rightcoord
            #4K
            #55hapfreqspec
            #6strainnum
            #7H1 (index 6)
            #8H2 (index 7)
            #9H12 (index 8)
            #10H2/H1 (index 9)
            #11H123 (index 10)
 
    f.close()

    return pos_list, G1_list

   

def rankGenes(ann_dict, chrom_dict):
    outfileName=save_prefix+'_g1_ranked.txt'
    G1_dict={}
    gene_G1s=[]
    for chrom in ann_dict:
        pos_list=chrom_dict[chrom][0]
        G1_list=chrom_dict[chrom][1]
        for gene in ann_dict[chrom]:
            G1s=[]
            start=ann_dict[chrom][gene][0]
            stop=ann_dict[chrom][gene][1]
            for i in range(len(pos_list)):
                if start<pos_list[i]<stop:
                    G1s.append(G1_list[i])
                    

            if len(G1s)>0:
                gene_G1=np.mean(G1s)
                gene_G1s.append(gene_G1)
                G1_dict[gene]=[gene_G1, chrom, start, stop]

    #print(sf2_dict['YAL054C'])

    with open(outfileName, 'w') as wf:
        wf.writelines('gene\taverageG1\tchrom\tstart\tstop\n')
        for k in sorted(G1_dict, key=G1_dict.get, reverse=True):
            #print(gene)
            gene=k
            avgG1=str(G1_dict[gene][0])
            chrom=str(G1_dict[gene][1])
            start=str(G1_dict[gene][2])
            stop=str(G1_dict[gene][3]) 
            wf.writelines(gene+'\t'+avgG1+'\t'+chrom+'\t'+start+'\t'+stop+'\n')

    return G1_dict,gene_G1s

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
    pos_list, G1_list=ParseH12(h12file)
    #get chromosome names
    prefix=h12file.split('.')[0]
    chrom=prefix.split('Spar_')[1]
    chrom_dict[chrom]=[pos_list, G1_list]

h12_dict, gene_G1s = rankGenes(ann_dict, chrom_dict)
#plotHistHstats(gene_G1s, 'G1')
