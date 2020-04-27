import os
import vcf
import sys


##USGAGE: python vcf2SweepFinderInput.py vcf_filename

##NOTE: NEED TO EDIT THIS TO WORK ON THE MERGED FILE INSTEAD, IT'S NOT QUITE RIGHT

folded = 0 #change for polarized data

def main():
        with open(outfile_name,"w") as wf:     
            vcf_reader = vcf.Reader(open(vcf_filename, 'r'))
            wf.writelines("position\tx\tn\tfolded\n")
            for record in vcf_reader:
                n=0 # allele count
                x=0 # derived allele count
                pos=record.POS
                ref=alleles_called=[record.alleles[0]]
                for call in record.samples:
                    if call.gt_type!=None: #ignoring uncalled sites 
                        n+=1
                        allele=call.gt_bases
                        if allele!=ref:
                                x+=1
                        
                
                if (folded==1 and x>0) or folded==0: #do not add invariant polarized sites
                        
                        wf.writelines(str(pos)+"\t"+str(x)+"\t"+str(n)+"\t"+str(folded)+"\n")
##                        print(str(pos)+"\t"+str(x)+"\t"+str(n)+"\t0\n")
##                        exit() #break for speediness
##                
                
                
if __name__ == '__main__':
        vcf_filename=sys.argv[1]
        outfile_name=vcf_filename.split('.')[0]+".alleleFreqFile"
        main()
