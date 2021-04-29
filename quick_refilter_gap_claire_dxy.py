import os
import pandas as pd
import numpy as np
from Bio import SeqIO
import csv
import sys

#USAGE:


#python3 dxy_calc.py alignment_directory output_file.csv
#alignment_directory contains one Scer population aligned with a Spar population, one file per gene

#command passed:
#python ~/scripts/quick_refilter_gap_claire_dxy.py ./ filtered_out_genes_dxy_gap0pt5.csv
#in /usr2/people/mabrams/data/claire_dxy_alignments_1011_EuroSpar/1011pops_EuroSpar_aligned/1WineEuropean


#filter out bad alignments
def separate_strains(file, gap_max=0.05, spar_strain_count_min=8, scer_strain_count_min=0): #.7*724
    
    spar_seqs, scer_seqs = [], []
    gap_ratios=[]

    for record in SeqIO.parse(file, 'fasta'):

        gap_ratio=(str(record.seq).count('-') + str(record.seq).count('N') + str(record.seq).count('n')) / len(record.seq)
        gap_ratios.append(gap_ratio)
        if gap_ratio > gap_max:

            continue

        if '_Sp_' in record.description:
            spar_seqs += [record.seq]

        elif 'ref(S288c)' not in record.description:

            scer_seqs += [record.seq]

    if len(scer_seqs) < scer_strain_count_min or len(spar_seqs) < spar_strain_count_min:

        return file[3:9],np.mean(gap_ratios),scer_seqs,spar_seqs

    return None,np.mean(gap_ratios),scer_seqs,spar_seqs


directory = sys.argv[1]
output = sys.argv[2]

if __name__ == '__main__':

    with open(output, 'w') as f:
        writer = csv.writer(f)
        writer.writerow(['directory', 'gene', 'avg_gap_ratio','len(spar_seqs)', 'len(scer_seqs)'])

        f.close()
                

    for file in os.listdir(directory):
        if "muscle_afa" not in file:
            continue
        gene = file.split('.')[0].split('_')[-1]

        filtered_out,avg_gap_ratio,scer_seqs,spar_seqs = separate_strains(directory+"/"+file)

        if filtered_out!=None:
            print(filtered_out,avg_gap_ratio)


            with open(output, 'a') as f:
                writer = csv.writer(f)

                writer.writerow([directory, gene, avg_gap_ratio,len(spar_seqs), len(scer_seqs)])

                f.close()


