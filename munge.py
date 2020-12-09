#Parameters

inFileName='geneFit_individual_replicates.txt'
outFileName='genFit_37_for_Tstat_MA.txt'
#outFileName='default'

#Functions

def munge(geneFit_individual_replicates_file,paired=True,outFileName='default'):
    '''input: geneFit_individual_replicates.txt file from RBSeq
       output: geneFit_temp_for_Tstat file.txt file for downstream analysis'''
    f=open(geneFit_individual_replicates_file)
    for line in f:
        row_data=line.strip().split('\t')
        if row_data[0]=='NearestGene': #if header row
            temp=row_data[3][:2] #get temp based on header content
            if outFileName=='default':
                outFileName='geneFit_'+temp+'_for_Tstat.txt'
            wf=open(outFileName,'w') #prepare new file with header
            wf_header='annotation\t'
            for biorep in row_data[3:]:
                br_fields=biorep.split('_')
                br_T1=br_fields[0]
                br_n1=br_fields[1]
                br_T2=br_fields[3]
                br_n2=br_fields[4]
                new_name='Biorep_'+br_n1+'_'+br_T1+'_vs_Biorep_'+br_n2+'_'+br_T2+'_fitness'
                wf_header+=new_name+'\t'
            wf_header+='gene\tallele\n'
            wf.writelines(wf_header)
        elif len(row_data)>3: #skip 'intergenic' line
            annotation=row_data[1]
            gene=annotation[2:]
            allele=annotation[:2]
            new_row=annotation+'\t'
            for biorep_data in row_data[3:]:
                new_row+=biorep_data+'\t'
            new_row+=gene+'\t'+allele+'\n'
            wf.writelines(new_row)
    wf.close()
    f.close()
    return None

##RUN##

munge(inFileName,outFileName=outFileName,paired=True)
                    
                
                

