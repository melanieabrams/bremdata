import sys

#Parameters
gff='/usr2/people/mabrams/Amended_Genomes/Z1/Z1.gff'
n=10 #top number to look at

#USAGE: python postrpocess_h12peaks.py *.h12peaks
#NOTE: program expects peakfiles to be named like Bergstrom2014_Spar_7.h12peaks

def ParseFromGFF(gfffile):
    '''
    Parses SGD features flat file
    Input: SGD_features.tab file
    Output: dict of {chrom:{gene:[start,stop]}}
    '''

    ann_dict={}
    
    f = open(gfffile)
    lines=[]
    for line in f:
        row_data=line.split("\t")
        chrom=row_data[0]
        start=int(row_data[3])
        stop=int(row_data[4])
        info=row_data[8].split(";")
        yName=info[0].split('=')[1]
        if chrom in ann_dict:
            ann_dict[chrom][yName]=[start,stop]
        else:
            ann_dict[chrom]={yName:[start,stop]}
    f.close()
    
    return ann_dict

def AnnotatePeakfile(peakfile, ann_dict):

    coding_lines=[]
    genes = []
    
    #get chromosome & outfile name from infile
    outfile_name=peakfile+'.annotated'
    prefix=peakfile.split('.')[0]
    chrom=prefix.split('_')[2]
    if len(chrom)==1:
        chrom='chr0'+chrom
    else:
        chrom='chr'+chrom
    print("chromosome: "+chrom)

    #write annotations to outfile
    headerline="1ctrcoord\t2leftcoord\t3rightcoord\t4K\t5hapfreqspec\t6strainnum\t7H1\t8H2\t9H12\10H2/H1\t11smedgepk\t12lgedgepk\t13gene\n"
    with open(outfile_name,"w") as wf:
        wf.writelines(headerline)
        f=open(peakfile)
        linecount=0
        for line in f:
            if linecount<n:
                row_data=line.split("\t")
                peakctr=int(row_data[0])
                #peakctr=6952#to test that it would correctly identify something in an orf chr01:6592 is in YAL062W
                #chrom='chr01'  #to test that it would correctly identify something in an orf: chr01:6592 is in YAL062W
                chrom_dict=ann_dict[chrom]
                #print("YAL062W" in chrom_dict) #to test that it would correctly identify something in an orf: chr01:6592 is in YAL062W
                annotation='noncoding'
                for gene in chrom_dict:
                    start=chrom_dict[gene][0]
                    stop=chrom_dict[gene][1]
                    if peakctr>start and peakctr<stop:
                        annotation=gene
                        genes.append(gene)
    ##                else:
    ##                    if gene=="YAL062W":
    ##                        print(start, type(start))
    ##                        print(stop, type(stop))
    ##                        print(peakctr, type(peakctr))
    ##                        print(peakctr>start)
    ##                        print(peakctr<stop)
    ##                        print(annotation)
    ##                        exit()
    ##            print(annotation)
    ##            exit()
                newline=line.strip('\n')+'\t'+annotation
                wf.writelines(newline+'\n')
                if annotation!='noncoding':
                    coding_lines.append(newline+'\t'+chrom+'\n')
                linecount+=1
        f.close()
    return coding_lines, genes

###START###
ann_dict=ParseFromGFF(gff)
peakfiles=sys.argv[1:]

all_coding_lines=[]
all_genes=[]
for peakfile in peakfiles:
    all_coding_lines+=AnnotatePeakfile(peakfile,ann_dict)[0]
    all_genes+=AnnotatePeakfile(peakfile,ann_dict)[1]

with open("all_coding_peaks.h12peaks.annotated",'w') as wf:
    for line in all_coding_lines:
        wf.writelines(line)
with open("genes_from_coding_peaks.txt",'w') as wf:
    for gene in all_genes:
        wf.writelines(gene+'\n')
