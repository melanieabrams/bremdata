import sys

#usage: python3 annotate_vcfs.py spar.gff file1.vcf file2.vcf...

def ParseFromGFF(gfffile):
    '''
    Parses SGD features flat file
    Input: SGD_features.tab file
    Output: dict of {gene:[chrom,start,stop]}
    '''

    ann_dict={}
    
    f = open(gfffile)
    lines=[]
    for line in f:
        row_data=line.split("\t")
        chrom=row_data[0]
        start=row_data[3]
        stop=row_data[4]
        col8=row_data[8].split(";")
        yName=col8[4][4:]
        ann_dict[yName]=[chrom,start,stop]
    f.close()
    
    return ann_dict

def getGenic(vcf_file,gff_dict):
    '''for each line, see if it is genic or intergenic based on gff'''
    #open vcf
    f=open(vcf_file)

    headerlines=[]
    varlines=[]
    for line in f:
        if line[0]=='#':
            headerlines.append(line)
        else:
            varlines.append(line)
    f.close()

    headerlines[-1]=headerlines[-1].strip("\n")+("\t"+"GENE"+"\n")
    
    #check if each variant call is genic
    filtered_lines=[]
    for line in varlines:
        chrom=line.split("\t")[0][2:]
        pos=int(line.split("\t")[1])
        for gene in gff_dict:
            if chrom==gff_dict[gene][0]:
                start=int(gff_dict[gene][1])
                stop=int(gff_dict[gene][2])
                if start <= pos <= stop:
                    genic_line=line.strip("\n")+("\t"+gene+"\n")
                    filtered_lines.append(genic_line)
    
    
    #write new file
    file_id=vcf_file.split('.')[0]
    wf=open(file_id+'.genic_annotated.vcf','w') #outfile
    for line in headerlines:
        wf.writelines(line)
    for line in filtered_lines:
        wf.writelines(line)

    wf.close()

    return
    

    
    

gff=sys.argv[1]
files=sys.argv[2:]

ann_dict=ParseFromGFF(gff)

for eachFile in files:
    getGenic(eachFile,ann_dict)

    


