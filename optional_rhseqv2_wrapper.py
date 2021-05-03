import os

## USAGE ##
# python optional_rhseqv2_wrapper
# optional wrapper for running the rhseq pipeline from poolCount to mwu analysis with one command
# generates, prints (optional), and runs (optional) a command for all the rhseq scripts with their current paramters
# call this wrapper from the folder where output is intended)

##PARAMETERS##
printCommand=True
runCommand=False

scriptDir='/usr2/people/mabrams/scripts/' #the location where the rhseqv2 scripts are saved
py3='python' #the call for python at the start of the script


poolCountFile='poolCount.txt'#name/path to the first infile (the poolCount from the Coradetti et al. RBseq pipeline RBseq_Count_Barcodes script v1.1.4 (no changes))


##RUN##
if __name__ == '__main__':
    command=''

    cmd_start=py3+' '+scriptDir
    outDir='./' 
    cmd_mid=' '+outDir+' '
    
    command=cmd_start+'convert_poolcount_to_fastqpooledreadsclean_rhseqv2.py'+cmd_mid+poolCountFile+';'
    command+=cmd_start+'2total_reads_and_normalize_rhseqv2.py'+cmd_mid+'./*fastq_pooled_reads_clean;'
    command+=cmd_start+'3remove_NC_and_plasmid_inserts_rhseqv2.py'+cmd_mid+'./*.normalized_pooled_reads;'
    command+=cmd_start+'4reformat_single_techrep_rhseqv2.py'+cmd_mid+'./*.normalized_pooled_reads_coding;'
    command+=cmd_start+'5combine_bioreps_rhseqv2.py'+cmd_mid+'*.normalized_averaged_techreps;'
    command+=cmd_start+'6filter_inserts_rhseqv2.py'+cmd_mid+'*.normalized_averaged_bioreps;'
    command+=cmd_start+'7fitness_ratios_rhseqv2.py'+cmd_mid+'*.filtered_inserts;'
    command+=cmd_start+'8organize_and_filter_genes_rhseqv2.py'+cmd_mid+'*.insert_ratios;'
    command+=cmd_start+'9mann_whitney_u_rhseqv2.py'+cmd_mid+'*.filtered_gene_inserts;'
    command+=cmd_start+'10effect_size_rhseqv2.py'+cmd_mid+'*.filtered_gene_inserts'+'*.mwu_test_results;'


    if printCommand==True:
        print(command)
    if runCommand==True:
        os.system(command)

    




