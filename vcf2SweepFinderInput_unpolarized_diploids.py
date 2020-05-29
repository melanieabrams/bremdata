import os
import vcf
import sys


##USGAGE: python vcf2SweepFinderInput.py vcf_filename

##PARAMETERS

#samples_to_exclude=['N-44','YPS138','UFRJ50816','IFO1804','Q31.4'] # for Bergstrom 2014 S paradoxus
        #the first four of those ('N-44','YPS138','UFRJ50816','IFO1804') are outside of the population of interest but in the source vcf
        #the last one (Q31.4) is there because it has a LOT of uncalled sites in the vcf

samples_to_exclude=[] # for 1011 genomes (all pop)

diploid_present=True # false for Bergstrom 2014 S paradoxus, true  for 1011 genomes S cerevisiae

folded = 1 #for unpolarized data, as hardcoded


def main():
        samples_tested=[]
        with open(outfile_name,"w") as wf:     
            vcf_reader = vcf.Reader(open(vcf_filename, 'r'))
            wf.writelines("position\tx\tn\tfolded\n")
            for record in vcf_reader:
                allele_dict={}                
                pos=record.POS
                n=0
                for call in record.samples:
                    sample=call.sample
                    if sample not in samples_to_exclude:
                            if sample not in samples_tested:
                                    samples_tested.append(sample)
                            gt=call.gt_bases
                            if call.gt_type!=None and gt!=None: #ignoring uncalled sites
                                    alleles=gt.split('/')
                                    for allele in alleles:
                                            if len(allele)==1: #ignore indels
                                                    n+=1 #add each allele to sample size
                                                    if allele in allele_dict:
                                                            allele_dict[allele]+=1
                                                    else:
                                                            allele_dict[allele]=1
                        
                #print(allele_dict)
                if allele_dict!={}: #ignore sites where no SNPs in this population
                        sample_alleles=allele_dict.values()
                        x=max(sample_alleles)
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
