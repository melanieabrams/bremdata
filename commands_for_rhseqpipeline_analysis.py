
###1

#first,make pooLCount file with sam's codebase
#next, edit the test temperature in convert_poolcount_to_fastqpooledreadsclean_2.py
#then, run format conversion script

#python ~/scripts/convert_poolcount_to_fastqpooledreadsclean_v2.py ./ poolCount.txt

###2
#~/bin/anaconda2/bin/python ~/carlygithub/rh-seq/2total_reads_and_normalize.py ./ ./*fastq_pooled_reads_clean

###3
#~/bin/anaconda2/bin/python ~/carlygithub/rh-seq/3remove_NC_and_plasmid_inserts.py ./ ./*.normalized_pooled_reads


###4
# ~/bin/anaconda2/bin/python ~/scripts/4one_tr.py ./ ./*.normalized_pooled_reads_coding




###5
##(had to modify this because py2.6 to 7)~/bin/anaconda2/bin/python ~/carlygithub/rh-seq/5combine_bio_reps.py ./ *.normalized_averaged_techreps
## modified the groupby statements in all further scripts
# ~/bin/anaconda2/bin/python ~/scripts/5combine_bio_reps_MBA.py ./ *.normalized_averaged_techreps


##for 38C:
##The mean coeffvar of  39  is  0.77385306228
##The mean coeffvar of  28  is  0.651366182018

###6-8 collapsed br

#~/bin/anaconda2/bin/python ~/scripts/6filter_inserts_uncollapsedMBA.py ./ *.normalized_averaged_bioreps



##also ran a version of 6 with coevar 2.0 instead of 1.5
#~/bin/anaconda2/bin/python ~/scripts/6filter_inserts_uncollapsedMBA_newsettings.py ./ *.normalized_averaged_bioreps


##548129 before filtering
##480055 after filtering for coeffvar
##168044 after filtering for coeffvar and read count

##~/bin/anaconda2/bin/python ~/scripts/7fitness_ratios_uncollapsedMBA.py ./ *.filtered_inserts
##~/bin/anaconda2/bin/python ~/scripts/8organize_and_filter_genes_ratiocv_uncollapsedMBA.py ./ 1.5_1.1_coeffvarANDreadcutoff_uncollapsedMBA.insert_ratios
##with 
##number_inserts_per_allele_needed_to_test =2.0
##cutoff_gene_cv_ratio = 20.0 #coefficient of variation of the ratio of 39/28 for an insert


##before filtering
##168044
##after filtering
##38341
##too few inserts: 2663
##too high cv: 12
##
##with


##9
#python /usr2/people/mabrams/scripts/9mann_whitney_u_uncollapsedMBA_py3_testgene.py ./ 1.5_1.1_coeffvarANDreadcutoff_uncollapsedMBA_ratiocvfilt_20.0_2.0.filtered_gene_inserts

##collapsed
#~/bin/anaconda2/bin/python ~/carlygithub/rh-seq/5combine_bio_reps.py ./ *.normalized_averaged_techreps
#~/bin/anaconda2/bin/python ~/carlygithub/rh-seq/6filter_inserts_bug_fixed_by_AF_py2pt7.py  ./ *.normalized_averaged_bioreps
#~/bin/anaconda2/bin/python ~/carlygithub/rh-seq/7fitness_ratios.py  ./ *.filtered_inserts
#~/bin/anaconda2/bin/python ~/carlygithub/rh-seq/8organize_and_filter_genes_py2pt7.py  ./ 1.5_1.1_coeffvarANDreadcutoff.insert_ratios 

#python ~/carlygithub/rh-seq/9mann_whitney_u_py3.py ./ 1.5_1.1_coeffvarANDreadcutoff_5.0.filtered_gene_inserts
