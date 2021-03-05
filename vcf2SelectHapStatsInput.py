import os
import vcf
import sys


##USGAGE: python vcf2SelectHapStatsInput.py vcf_filename


##PARAMETERS
#samples_to_exclude=['N-44','YPS138','UFRJ50816','IFO1804','Q31.4'] # for Bergstrom 2014 S paradoxus
        #the first four of those ('N-44','YPS138','UFRJ50816','IFO1804') are outside of the population of interest but in the source vcf
        #the last one (Q31.4) is there because it has a LOT of uncalled sites in the vcf

samples_to_exclude=[] # for 1011 genomes (all pop)

heterozygous_present=True # false for Bergstrom 2014 S paradoxus, true  for 1011 genomes S cerevisiae



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
                                if heterozygous_present==True:
                                        bases=allele.split('/')
                                        if bases[0]!=bases[1]:
                                                allele='.' # represent heterozygous GT as dots
                                        else:
                                                allele=bases[0] #if homozygous, represent GT as the base
                                if len(allele)>1:
                                        allele = 'INDEL'
                            else:
                                allele='None'
                                if sample in missing_data:
                                        missing_data[sample]+=','+str(pos)
                                else:
                                        missing_data[sample]=str(pos)
                            pos_line+=","+str(allele)
                if 'INDEL' not in pos_line and 'None' not in pos_line:
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
