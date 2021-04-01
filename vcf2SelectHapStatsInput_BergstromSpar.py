import os
import vcf
import sys


##USGAGE: python vcf2SelectHapStatsInput.py vcf_filename
##NOTE: this script assumes standard VCF formatting for the vcf file.  Unlike for gVCF, it handles missing data as reference.


##PARAMETERS


samples_to_exclude=[] # for for Bergstrom et al., S.paradoxus all populations

heterozygous_present=False # True for Bergstrom et al., S.paradoxus



##RUN


class MyStdErr(object):

    stderr = sys.stderr
    waserr = False

    def write(self, text):
        self.waserr = True
        self.stderr.write(text)




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
                            allele=record.REF #assume reference till variant for this vcf format
                            try:
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
##                            else:
##                                allele='None'
##                                if sample in missing_data:
##                                        missing_data[sample]+=','+str(pos)
##                                else:
##                                        missing_data[sample]=str(pos)
                            except:
                                    None
                            pos_line+=","+str(allele)
                if 'INDEL' not in pos_line:# and 'None' not in pos_line:
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
        sys.stderr = MyStdErr()
        vcf_filename=sys.argv[1]
        outfile_name=vcf_filename.split('.')[0]+".shsinput"
        errfile_name=vcf_filename.split('.')[0]+".missing_call_log"
        main()
