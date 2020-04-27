import os
import vcf
import sys


##USGAGE: python vcf2SelectHapStatsInput.py vcf_filename

##NOTE: NEED TO EDIT THIS TO WORK ON THE MERGED FILE INSTEAD, IT'S NOT QUITE RIGHT



def main():
        with open(outfile_name,"w") as wf:
            vcf_reader = vcf.Reader(open(vcf_filename, 'r'))
            for record in vcf_reader:
                pos_line=""
                pos=record.POS
                pos_line+=str(pos)
                for call in record.samples:
                    if call.gt_type!=None:
                        allele=call.gt_bases
                    else:
                        allele='.'
                    pos_line+=","+str(allele)
                wf.writelines(pos_line+"\n")
##                print(pos_line)
##                exit() #break for speediness
                
                
                
if __name__ == '__main__':
        vcf_filename=sys.argv[1]
        outfile_name=vcf_filename.split('.')[0]+".shsinput"
        main()
