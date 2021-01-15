import sys

from resample_g1 import * #both in the ~/scripts/ folder of /usr2/people/mabrams on thar.  Resample g1 and g12 are the same script, pulling the same column of their respective rank files, g2g1 pulls from the second column.


#same as automate-g1_test scripts but use the ranked file with g12 as the first column after gene and g2g1 as second

##Parameters##
essential_gene_path='/usr2/people/mabrams/Amended_Genomes/essential.csv'
gff='/usr2/people/mabrams/Amended_Genomes/S288C/saccharomyces_cerevisiae_R64-1-1_20110208_withoutFasta.gff'
printLots=True #set to off for less output text 

ranked_g12_paths={
    '1Wine':'/usr2/people/mabrams/data/g1_5top/1Wine/merged_1WineEuropean_ranked.txt',
    'M3': '/usr2/people/mabrams/data/g1_5top/M3/merged_M3_MosaicRegion3_ranked.txt',
    '8Mixed': '/usr2/people/mabrams/data/g1_5top/8Mixed/merged_8MixedOrigin_ranked.txt',
    '25Sake': '/usr2/people/mabrams/data/g1_5top/25Sake/merged_25Sake_ranked.txt',
    '3Brazil': '/usr2/people/mabrams/data/g1_5top/3Brazil/merged_3BrazilianBioethanol_ranked.txt'}

##goi_lists={
##    'published analysis of Weiss et al. 2018 RH-Seq, padj<0.05':[
##        'YLR397C','YGR098C', 'YMR168C','YKR054C',
##        'YDR180W','YHR023W','YCR042C','YNL172W'],
##    'paired reanalysis of Weiss et al. 2018 RH-Seq':[
##        'YLR397C','YGR098C', 'YMR168C','YKR054C',
##        'YDR180W','YHR023W','YGR198W','YHR166C',
##        'YCR042C','YPL174C'],
##    'unpaired reanalysis of Weiss et al. 2018 RH-Seq, median, padj<0.1':[
##        'YLR397C','YGR098C','YKR054C','YDR180W',
##        'YCR042C','YMR168C','YHR166C','YGR198W',
##        'YHR023W','YNL172W','YHR204W','YJR135C',
##        'YER090W','YBL066C','YJL047C','YKR068C',
##        'YNL238W','YGR246C'],
##    'unpaired reanalysis of Weiss et al. 2018 RH-Seq, median, padj<0.05 or 0.01':[
##        'YLR397C','YGR098C','YKR054C','YDR180W',
##        'YCR042C','YMR168C','YHR166C','YGR198W',
##        'YHR023W','YNL172W','YHR204W','YJR135C',
##        'YER090W','YBL066C'],
##    'unpaired reanalysis of Weiss et al. 2018 RH-Seq, median, p_adj<0.001':[
##        'YLR397C','YGR098C','YKR054C','YDR180W',
##        'YCR042C','YMR168C','YHR166C','YGR198W'],
##    'unpaired, reanalysis of Weiss et al. 2018 RH-Seq, mean, p_adj<0.1 - excluding YPR164W b/c no g1':[
##        'YLR397C','YGR098C','YMR168C','YKR054C',
##        'YHR023W','YDR180W','YPL174C','YCR042C',
##        'YMR016C','YJR135C','YJL025W','YDR443C',
##        'YKL134C','YGR140W','YHR166C'],
##    'unpaired, reanalysis of Weiss et al. 2018 RH-Seq, mean, p_adj<0.05 - excluding YPR164W b/c no g1':[
##        'YLR397C','YGR098C','YMR168C','YKR054C',
##        'YHR023W','YDR180W','YPL174C','YCR042C',
##        'YMR016C','YJR135C','YJL025W','YDR443C',
##        'YKL134C'],
##    'unpaired, reanalysis of Weiss et al. 2018 RH-Seq, mean, p_adj<0.01 - excluding YPR164W b/c no g1':[
##        'YLR397C','YGR098C','YMR168C','YKR054C',
##        'YHR023W','YDR180W','YPL174C','YCR042C'],
##    'unpaired, reanalysis of Weiss et al. 2018 RH-Seq, mean, p_adj<0.001 - excluding YPR164W b/c no g1':[
##        'YLR397C','YGR098C','YMR168C','YKR054C',
##        'YHR023W','YDR180W'],
##    }

goi_lists={

    'unpaired, reanalysis of Weiss et al. 2018 RH-Seq, mean, p_adj<0.05 - excluding YPR164W b/c no g1':[
        'YLR397C','YGR098C','YMR168C','YKR054C',
        'YHR023W','YDR180W','YPL174C','YCR042C',
        'YMR016C','YJR135C','YJL025W','YDR443C',
        'YKL134C'],
    }





##Functions##

def getResample(pop,goi_list,gff_genes):
    if printLots==True:
        print('\n\npopulation:\t'+pop)
    ranked_genes_G12=ranked_g12_paths[pop]
    tested_goi=test_goi_in_gff(goi_list,gff_genes,printGOI=printLots)
    essential_genes=ParseEssential(essential_gene_path)
    G12_dict=ParseRankedGenes(ranked_genes_G12,goi_list,printGeneG1=printLots)
    rv=resample(G12_dict,essential_genes,goi_list,printResampling=printLots) #resample that given G12 and don't print longform (gene by gene, num essential, etc)
    if rv<0.05:
        print(pop,str(rv)+'*')
        return rv, 's'
    else:
        print(pop, rv)
        return rv, 'ns'

def resampleAllCombos(ranked_g12_paths,goi_lists,gff_genes):
    for goi_list_desc in goi_lists:
        significant=0
        print('==='+goi_list_desc+'===')
        goi_list=goi_lists[goi_list_desc]
        for pop in ranked_g12_paths:
            p,sig=getResample(pop,goi_list,gff_genes)
            if sig=='s':
                significant+=1
        print('significant (rv<0.05): ', significant)
            


###RUN###
if __name__ == "__main__":
    
    gff_genes=ParseFromGFF(gff)
    resampleAllCombos(ranked_g12_paths,goi_lists,gff_genes)
    
