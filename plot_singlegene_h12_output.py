import sys
import numpy as np
import matplotlib.pyplot as plt

'''
rank selecthapstats h12 outputs by gene
assuming file name in form: 
merged_6AfricanBeer_chromosome1.h12_h2h1

'''
#Parameters
#gff='/usr2/people/mabrams/Amended_Genomes/D1373/DBVPG1373.gff'
gff='/usr2/people/mabrams/Amended_Genomes/S288C/saccharomyces_cerevisiae_R64-1-1_20110208_withoutFasta.gff'
roman_numerals_in_gff=True

goi=['YLR397C','YGR098C', 'YMR168C','YKR054C','YDR180W',
     'YHR023W','YGR198W','YHR166C','YCR042C','YPL174C']

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

   

def getGeneStats(ann_dict, chrom_dict):
    '''from rank_genes function in ranking output script
    returns a dictionary of gens with their h12s, h12 ratios, chrom, start, and stop'''
    H12_dict={}
    for chrom in ann_dict:
        pos_list=chrom_dict[chrom][0]
        H12_list=chrom_dict[chrom][1]
        ratio_H2H1_list=chrom_dict[chrom][2]
        for gene in ann_dict[chrom]:
            H12s=[]
            H2H1_ratios=[]
            positions=[]
            start=ann_dict[chrom][gene][0]
            stop=ann_dict[chrom][gene][1]
            for i in range(len(pos_list)):
                if start<pos_list[i]<stop:
                    positions.append(int(pos_list[i]))
                    H12s.append(H12_list[i])
                    H2H1_ratios.append(ratio_H2H1_list[i])            
            H12_dict[gene]=[H12s, H2H1_ratios, positions, chrom, start, stop]

    return H12_dict

def plotStatOverGene(gene, statdict, stat='h12'):
    if stat=='h12':
        y=np.array(statdict[gene][0])
    else:
        y=np.array(statdict[gene][1])
    x=np.array(statdict[gene][2])
    chrom=statdict[gene][3]
    start=statdict[gene][4]
    stop=statdict[gene][5]
    
    fig=plt.figure(figsize=(20,15))
    pl = plt.subplot(111)

    pl.scatter(x,y)
    
    pl.set_xlabel('pos',fontsize=50)
    pl.set_ylabel(stat,fontsize=50)

    plt.title(prefix+' '+stat+': '+gene+'('+chrom+':'+str(start)+'-'+str(stop)+')')
    plt.savefig(save_prefix+stat+'_'+gene+'.pdf', bbox_inches='tight', format='pdf', dpi=1000)

    plt.clf()
    return None



##RUN###

ann_dict=ParseFromGFF(gff)
chrom_dict={}

#build chrom dict
print('reading infiles')
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

print('getting gene stats')
h12_dict= getGeneStats(ann_dict, chrom_dict)

print('plotting')
for gene in goi:
    plotStatOverGene(gene,h12_dict,'h12')
    
