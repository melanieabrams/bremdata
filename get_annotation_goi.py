import urllib.parse

##PARAMETERS##
gff='saccharomyces_cerevisiae_R64-1-1_20110208_withoutFasta.gff'


goi=['YGR198W','YLR397C','YBR136W','YGR140W','YNL049C','YKL017C','YGL082W',
     'YMR125W','YDR508C','YMR207C','YDR375C','YDR180W','YKL197C','YDR318W',
     'YMR094W','YOR326W','YBR081C','YPR049C','YIL152W','YER151C','YJR107W',
     'YAL026C','YDR456W','YLR141W','YPL268W','YDR235W']

outFile='barseq_goi_sgd_notes.tsv'

##FUNCTIONS##
def ParseFromGFF(gff):
    '''
    Parses gff
    Output: dict of {yName:{geneName}}
    '''
    ann_dict={}
    
    f = open(gff)
    lines=[]
    for line in f:
        if line[0]!='#': #skip header rows
            row_data=line.split("\t")
            if row_data[2]=='gene':
                info=row_data[8].split(";")
                yName=info[0].split('=')[1]
                ann_dict[yName]=[yName,'no desc']
                for i in info:
                    if i.startswith('gene'):
                        geneName=i.split('=')[1]
                        ann_dict[yName][0]=geneName
                    if i.startswith('Note'):
                     geneNote=urllib.parse.unquote(i.split('=')[1])
                     ann_dict[yName][1]=geneNote
                        
    f.close()
    return ann_dict


##RUN##
gff_dict=ParseFromGFF(gff)
with open(outFile,'w') as wf:
    wf.writelines('BarSeqHit\tgeneName\tSGD_note\n')
    for gene in goi:
    #    None
        print(gff_dict[gene])
        wf.writelines(gene+'\t'+gff_dict[gene][0]+'\t'+gff_dict[gene][1]+'\n')
    

    
