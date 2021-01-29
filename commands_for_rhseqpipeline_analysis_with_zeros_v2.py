
###1

#first,make pooLCount file with sam's codebase
#next, edit the test temperature in convert_poolcount_to_fastqpooledreadsclean_2.py
#then, run format conversion script

#then, RHseq pipeline
#python ~/scripts/convert_poolcount_to_fastqpooledreadsclean_37C.py ./ poolCount.txt ; ~/bin/anaconda2/bin/python ~/carlygithub/rh-seq/2total_reads_and_normalize.py ./ ./*fastq_pooled_reads_clean; ~/bin/anaconda2/bin/python ~/carlygithub/rh-seq/3remove_NC_and_plasmid_inserts.py ./ ./*.normalized_pooled_reads; ~/bin/anaconda2/bin/python ~/scripts/4one_tr.py ./ ./*.normalized_pooled_reads_coding; ~/bin/anaconda2/bin/python ~/scripts/5combine_bio_reps_MBA_nonzeros.py ./ *.normalized_averaged_techreps ; ~/bin/anaconda2/bin/python ~/scripts/6filter_inserts_uncollapsedMBA_debug.py ./ *.normalized_averaged_bioreps; ~/bin/anaconda2/bin/python ~/scripts/7fitness_ratios_uncollapsedMBA_nozeros_unpaired.py ./ *.filtered_inserts; ~/bin/anaconda2/bin/python ~/scripts/8organize_and_filter_genes_ratiocv_uncollapsedMBA_debug_nozeros.py ./ *.insert_ratios; python /usr2/people/mabrams/scripts/9mann_whitney_u_uncollapsedMBA_py3_testgene_nonzeros.py ./ *.filtered_gene_inserts
