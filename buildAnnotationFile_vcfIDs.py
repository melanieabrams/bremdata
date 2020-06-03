import os
import vcf
import sys


##USGAGE: python builtAnnotationFile_vcfIDs.py vcf_filename


##PARAMETERS


##RUN


def main():
       
        with open(outfile_name,"w") as wf:
            vcf_reader = vcf.Reader(open(vcf_filename, 'r'))
            for record in vcf_reader:
                chrom=str(record.CHROM)
                pos=str(record.POS)
                new_ID='chr'+chrom+'_'+pos
                wf.writelines(chrom+'\t'+pos+'\t'+new_ID+'\n')
                
        
                
if __name__ == '__main__':
        vcf_filename=sys.argv[1]
        outfile_name=vcf_filename.split('.')[0]+"_IDs.txt"
        main()
