import os
import vcf
import sys


##USGAGE: python vcf2SelectHapStatsInput.py vcf_filename

##NOTE: NEED TO EDIT THIS TO WORK ON THE MERGED FILE INSTEAD, IT'S NOT QUITE RIGHT

##PARAMETERS
samples_to_exclude=['N-44','YPS138','UFRJ50816','IFO1804','Q31.4']
#the first four of those ('N-44','YPS138','UFRJ50816','IFO1804') are outside of the population of interest but in the source vcf
#the last one (Q31.4) is there because it has a LOT of uncalled sites in the vcf

##RUN


def main():
        missing_data={} #samples with uncalled variant positions
        samples_tested=""
        num_samples_tested=0

        ##
        with open(outfile_name,"w") as wf:
            vcf_reader = vcf.Reader(open(vcf_filename, 'r'))
            for record in vcf_reader:
                pos_line=""
                pos=record.POS
                pos_line+=str(pos)
                for call in record.samples:
                    sample=call.sample
                    if sample not in samples_to_exclude:
                            if sample not in samples_tested:
                                    samples_tested+=sample+", "
                                    num_samples_tested+=1
                            if call.gt_type!=None and call.gt_bases!=None:
                                allele=call.gt_bases
                            else:
                                allele='.'
                                if sample in missing_data:
                                        missing_data[sample]+=','+str(pos)
                                else:
                                        missing_data[sample]=str(pos)
                            pos_line+=","+str(allele)
                if '.' not in pos_line:
                        wf.writelines(pos_line+"\n") #omit lines with missing data
##                print(pos_line)
##                exit() #break for speediness

        with open(errfile_name,"w") as wf:
                wf.writelines("the following samples were excluded")
                wf.writelines("%s," % sample for sample in samples_to_exclude)
                wf.writelines("shs input generated for the following"+str(num_samples_tested)+"samples\n")
                wf.writelines(samples_tested+"\n")
                wf.writelines("sample\tuncalled_pos\n")
                for sample in missing_data:
                        wf.writelines(sample+"\t"+missing_data[sample]+"\n")
                
        
                
if __name__ == '__main__':
        vcf_filename=sys.argv[1]
        outfile_name=vcf_filename.split('.')[0]+".shsinput"
        errfile_name=vcf_filename.split('.')[0]+".missing_call_log"
        main()
