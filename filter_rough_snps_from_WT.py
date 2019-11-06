import os
import vcf
import argparse


def main():
        wt_reader = vcf.Reader(open(args.w, 'r'))
        print("hello")
        vcf_reader = vcf.Reader(open(args.i, 'r'))
        vcf_writer = vcf.Writer(open(args.o, 'w'), vcf_reader)

        tolerance=args.t
        
        wt_dict={}
        for record in wt_reader:
            chrom=record.CHROM
            pos=record.POS
            if chrom in wt_dict.keys():
                wt_dict[chrom].append(pos)
            else:
                wt_dict[chrom]=[pos]
                
        for record in vcf_reader:
            chrom=record.CHROM
            pos=record.POS
            if chrom not in wt_dict.keys():
                vcf_writer.write_record(record)
            else:
                too_close=False
                for i in wt_dict[chrom]:
                    if abs(int(i)-int(pos)) <= int(tolerance):
                        too_close=True
                if too_close == False:
                    vcf_writer.write_record(record)
                

if __name__ == '__main__':
        # argument parser
        parser = argparse.ArgumentParser(description='Filter out variants found in WT')
        # input/output

        parser.add_argument('-w', required = True, help = 'wt vcf')
        parser.add_argument('-i', required = True, help = 'input vcf')
        parser.add_argument('-o', required = True, help = 'output vcf')
        parser.add_argument('-t', required = True, help = 'tolerance (bp)')
        
        args = parser.parse_args()
        main()


