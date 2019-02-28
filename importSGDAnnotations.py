

def ParseFromFile(sgdflatfile):
    '''
    Parses SGD features flat file
    Input: SGD_features.tab file
    Output: 
    '''

    descr_dict={}
    
    f = open(sgdflatfile)
    for line in f:
        row_data =line.split("\t")
        yName=row_data[3]
        descr=row_data[15]
        descr_dict[yName]=descr
 #   print descr_dict["YGR098C"]
    return descr_dict
    
    
    

##RUN

SGD_flat=ParseFromFile("/Users/Melanie/Desktop/Melanie/NGS/SGD_features.tab")
#gene_annotation=get_annotation(gene)
#print gene_annotation
