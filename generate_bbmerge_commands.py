##BBMERGE
##cmds=''
##for i in range(1,11):
##    lib=str(i+26)
##    snum=str(i)
##    cmd='~/bin/bbmap/bbmerge.sh in1=RBMA0'+lib+'A_S'+snum+'_L001_R1_001.fastq in2=RBMA0'+lib+'A_S'+snum+'_L001_R2_001.fastq out='+lib+'_merged.fq outu='+lib+'_unmerged.fq ihist='+lib+'ihist.txt; '
##    #print(cmd)
##    cmds+=cmd
##print(cmds)

###MAP AND POOL BLAT
##cmds=''
##for i in range(27,37):
##    merged=str(i)+'_merged.fq '
##    cmd='/usr/bin/python ~/carlygithub/rh-seq/map_and_pool_BLAT_2-1-17.py '+merged+'./ ; '
##    #print(cmd)
##    cmds+=cmd
##print(cmds)

#####1
##cmds=''
##for i in range(27,37):
##    fname=str(i)+'_merged.fq_pooled_reads '
##    cmd='/usr/bin/python ~/carlygithub/rh-seq/1clean_JR_pooled_output.py ~/carlygithub/rh-seq/YS2+CBS433+plasmid_clean ./'+fname+'; '
##    #print(cmd)
##    cmds+=cmd
##print(cmds)



##manually v'd library names to temp/br name (28A1.fq_clean etc), 28 and 39 even though 28 and 37, reran 2-4
##NOTE: Discovered that step 2-3 introduces a bug (naming reps like 280.0 not 28E1) if the tech rep is E.  make them H instead.


###2
#~/bin/anaconda2/bin/python ~/carlygithub/rh-seq/2total_reads_and_normalize.py ./ ./*.fq_clean

###3
#~/bin/anaconda2/bin/python ~/carlygithub/rh-seq/3remove_NC_and_plasmid_inserts.py ./ ./*.normalized_pooled_reads


###4
# ~/bin/anaconda2/bin/python ~/scripts/4one_tr.py ./ ./*.normalized_pooled_reads_coding






###5
##(had to modify this because py2.6 to 7)~/bin/anaconda2/bin/python ~/carlygithub/rh-seq/5combine_bio_reps.py ./ *.normalized_averaged_techreps
## modified the groupby statements in all further scripts
# ~/bin/anaconda2/bin/python ~/scripts/5combine_bio_reps_MBA.py ./ *.normalized_averaged_techreps


###6-8 collapsed br

#~/bin/anaconda2/bin/python ~/scripts/6filter_inserts_uncollapsedMBA.py ./ *.normalized_averaged_bioreps
##also ran a version of 6 with coevar 2.0 instead of 1.5
#~/bin/anaconda2/bin/python ~/scripts/6filter_inserts_uncollapsedMBA_newsettings.py ./ *.normalized_averaged_bioreps

##~/bin/anaconda2/bin/python ~/scripts/7fitness_ratios_uncollapsedMBA.py ./ *.filtered_inserts
##~/bin/anaconda2/bin/python ~/scripts/8organize_and_filter_genes_ratiocv_uncollapsedMBA.py ./ 1.5_1.1_coeffvarANDreadcutoff_uncollapsedMBA.insert_ratios
##with 
##number_inserts_per_allele_needed_to_test =2.0
##cutoff_gene_cv_ratio = 20.0 #coefficient of variation of the ratio of 39/28 for an insert
##
##with


##9
#python /usr2/people/mabrams/scripts/9mann_whitney_u_uncollapsedMBA_py3_testgene.py ./ 1.5_1.1_coeffvarANDreadcutoff_uncollapsedMBA_ratiocvfilt_20.0_2.0.filtered_gene_inserts

##uncollapsed
#~/bin/anaconda2/bin/python ~/carlygithub/rh-seq/5combine_bio_reps.py ./ *.normalized_averaged_techreps
#~/bin/anaconda2/bin/python ~/carlygithub/rh-seq/6filter_inserts_bug_fixed_by_AF_py2pt7.py  ./ *.normalized_averaged_bioreps
#~/bin/anaconda2/bin/python ~/carlygithub/rh-seq/7fitness_ratios.py  ./ *.filtered_inserts
#~/bin/anaconda2/bin/python ~/carlygithub/rh-seq/8organize_and_filter_genes_py2pt7.py  ./ 1.5_1.1_coeffvarANDreadcutoff.insert_ratios 

#python ~/carlygithub/rh-seq/9mann_whitney_u_py3.py ./ 1.5_1.1_coeffvarANDreadcutoff_5.0.filtered_gene_inserts
