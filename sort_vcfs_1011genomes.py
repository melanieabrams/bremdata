import os
poptable_file="1011_pops_poptable.csv"
vcf_dir="~/popgen/not_European_wine/"

def parse_poptable(poptable_file):
    '''parse a csv file with
        -column of Standardized name as A
        -colum of Clades as B
        data from 1011 yeast genomes paper
        remove non-alphanumeric characters'''
    pop_dict={}
    with open(poptable_file) as f:
        next(f)
        for line in f:
            row_data=line.split(",") 
            strain=row_data[0]
            pop=row_data[1].split("(")[0].strip()
            pop=pop.replace(" ","")
            pop=pop.replace(".","")
            pop=pop.replace("/","")
            
            if pop=="":
                pop="no_pop_assigned"
            if pop in pop_dict:
                pop_dict[pop].append(strain)
            else:
                pop_dict[pop]=[strain]
    return pop_dict

def sort_pops(pop_dict):
    for pop in pop_dict:
        if os.path.isdir("pop")!=True: #make directory if it doesn't exist
            cmd = 'mkdir '+vcf_dir+pop
            #print(cmd)
            os.system(cmd)
        for strain in pop_dict[pop]:
            cmd = 'mv '+vcf_dir+strain+'.vcf* '+vcf_dir+pop+'/'
            #print(cmd)
            #exit()
            os.system(cmd)
            
pop_dict=parse_poptable(poptable_file)
sort_pops(pop_dict)



