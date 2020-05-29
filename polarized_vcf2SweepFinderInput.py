import os
import vcf
import sys


##USGAGE: python vcf2SweepFinderInput.py vcf_filename

##PARAMETERS

#samples_to_exclude=['N-44','YPS138','UFRJ50816','IFO1804','Q31.4'] # for Bergstrom 2014 S paradoxus
        #the first four of those ('N-44','YPS138','UFRJ50816','IFO1804') are outside of the population of interest but in the source vcf
        #the last one (Q31.4) is there because it has a LOT of uncalled sites in the vcf

samples_to_exclude=[] # for 1011 genomes (all pop)

heterozygous_present=True # false for Bergstrom 2014 S paradoxus, true  for 1011 genomes S cerevisiae

folded = 0 #change for polarized data

def main():
        if heterozygous_present==True:
                allele_len=3
        else:
                allele_len=1
        samples_tested=[]
        with open(outfile_name,"w") as wf:     
            vcf_reader = vcf.Reader(open(vcf_filename, 'r'))
            wf.writelines("position\tx\tn\tfolded\n")
            for record in vcf_reader:
                n=0 # allele count
                x=0 # derived allele count
                pos=record.POS
                ref=record.alleles[0]
                if heterozygous_present==True:
                        ref=ref+'/'+ref
                for call in record.samples:
                    sample=call.sample
                    if sample not in samples_to_exclude:
                            if sample not in samples_tested:
                                    samples_tested.append(sample)
                            allele=call.gt_bases
                            if call.gt_type!=None and allele!=None: #ignoring uncalled sites
                                     if len(allele)==allele_len and len(ref)==allele_len: #ignore indels
                                        n+=1
                                        if allele!=ref:
                                                x+=1
        ##                                        print(allele,ref)
        ##                                        break
                
                if n>0: # don't take sites where there were no SNPs 
                        if (folded==0 and x>0) or folded==1: #do not add invariant polarized sites or indels
                        
                                wf.writelines(str(pos)+"\t"+str(x)+"\t"+str(n)+"\t"+str(folded)+"\n")
        ##                        print(str(pos)+"\t"+str(x)+"\t"+str(n)+"\t0\n")
        ##                        exit() #break for speediness
                        


        with open(errfile_name,"w") as wf:
                wf.writelines(str(len(samples_tested))+" samples tested: \n")
                wf.writelines(', '.join(samples_tested))
                wf.writelines("\n\nsamples excluded\n")
                wf.writelines(', '.join(samples_to_exclude))
        
        
##                
                
                
if __name__ == '__main__':
        vcf_filename=sys.argv[1]
        outfile_name=vcf_filename.split('.')[0]+".alleleFreqFile"
        errfile_name=vcf_filename.split('.')[0]+".log_alleleFreqFile"
        main()
