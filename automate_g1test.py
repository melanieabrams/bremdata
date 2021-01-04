import sys

from resample_g1 import * #both in the ~/scripts/ folder of /usr2/people/mabrams on thar



##Parameters##
essential_gene_path='/usr2/people/mabrams/Amended_Genomes/essential.csv'
gff='/usr2/people/mabrams/Amended_Genomes/S288C/saccharomyces_cerevisiae_R64-1-1_20110208_withoutFasta.gff'


ranked_g1_paths={
    '1Wine':'/usr2/people/mabrams/data/g1_5top/1Wine/merged_1WineEuropean_g1_ranked.txt',
    'M3': '/usr2/people/mabrams/data/g1_5top/M3/merged_M3_MosaicRegion3_g1_ranked.txt',
    '8Mixed': '/usr2/people/mabrams/data/g1_5top/8Mixed/merged_8MixedOrigin_g1_ranked.txt',
    '25Sake': '/usr2/people/mabrams/data/g1_5top/25Sake/merged_25Sake_g1_ranked.txt',
    '3Brazil': '/usr2/people/mabrams/data/g1_5top/3Brazil/merged_3BrazilianBioethanol_g1_ranked.txt'}

goi_lists={
    'unpaired analysis of Weiss et al. 2018 RH-Seq, padj<0.05':['YLR397C','YGR098C', 'YMR168C','YKR054C','YDR180W','YHR023W','YGR198W','YHR166C','YCR042C','YNL172W','YHR204W','YJR135C','YER090W','YBL066C'],
    'unpaired analysis of Weiss et al. 2018 RH-Seq, top 10 genes':['YLR397C','YGR098C', 'YMR168C','YKR054C','YDR180W','YHR023W','YGR198W','YHR166C','YCR042C','YNL172W'],
    'unpaired analysis of Weiss et al. 2018 RH-Seq, p_adj<0.001':['YLR397C','YGR098C', 'YMR168C','YKR054C','YDR180W','YHR023W','YGR198W','YHR166C','YCR042C'],
    'unpaired, mean - YPR164W not in ranked genes':['YLR397C','YGR098C', 'YMR168C','YKR054C','YDR180W','YHR023W','YCR042C','YPL174C','YCR042C'],
    'Weiss analysis goi':['YLR397C','YGR098C', 'YMR168C','YKR054C','YDR180W','YHR023W','YCR042C','YNL172W']
    }

##Functions##

def getResample(pop,goi_list,gff_genes):
    ranked_genes_G1=ranked_g1_paths[pop]
    tested_goi=test_goi_in_gff(goi_list,gff_genes,printGOI=False)
    essential_genes=ParseEssential(essential_gene_path)
    G1_dict=ParseRankedGenes(ranked_genes_G1,printGeneG1=False)
    rv=resample(G1_dict,essential_genes,printResampling=False) #resample that given G1 and don't print longform (gene by gene, num essential, etc)
    if rv<0.05:
        print(pop,str(rv)+'*')
        return rv, 's'
    else:
        print(pop, rv)
        return rv, 'ns'

def resampleAllCombos(ranked_g1_paths,goi_lists,gff_genes):
    for goi_list_desc in goi_lists:
        significant=0
        print('==='+goi_list_desc+'===')
        goi_list=goi_lists[goi_list_desc]
        for pop in ranked_g1_paths:
            p,sig=getResample(pop,goi_list,gff_genes)
            if sig=='s':
                significant+=1
        print('significant', significant)
            


###RUN###
if __name__ == "__main__":
    
    gff_genes=ParseFromGFF(gff)
    resampleAllCombos(ranked_g1_paths,goi_lists,gff_genes)
    
