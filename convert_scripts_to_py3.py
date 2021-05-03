
newname_dict={'2total_reads_and_normalize.py':'2total_reads_and_normalize_rhseqv2.py',
             '3remove_NC_and_plasmid_inserts.py': '3remove_NC_and_plasmid_inserts_rhseqv2.py',
             '4one_tr.py':'4reformat_single_techrep_rhseqv2.py',
             '5combine_bio_reps_MBA_nonzeros.py':'5combine_bioreps_rhseqv2.py',
             '6filter_inserts_uncollapsedMBA_debug.py':'6filter_inserts_rhseqv2.py',
             '7fitness_ratios_uncollapsedMBA_nozeros_unpaired.py':'7fitness_ratios_rhseqv2.py',
             '8organize_and_filter_genes_ratiocv_uncollapsedMBA_debug_nozeros.py':'8organize_and_filter_genes_rhseqv2.py',
             '9mann_whitney_u_uncollapsedMBA_py3_testgene_nonzeros.py':'9mann_whitney_u_rhseqv2.py'}

#newname_dict={'2total_reads_and_normalize.py':'2total_reads_and_normalize_rhseqv2.py'}


replace_dict={'28':'ctrl','39':'exptl','.iteritems(':'.items(','print ':'print('}



files=newname_dict.keys()
for file in files:
    f=open(file)
    with open(newname_dict[file],'w') as wf:
        for line in f:
            if 'print ' in line:
                line=line[:-1]
                line+=')\n'
            for item in replace_dict:
                line=line.replace(item,replace_dict[item])
            wf.writelines(line)
    f.close()
