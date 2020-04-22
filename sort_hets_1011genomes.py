import os

##INPUTS
poptable_file="1011_pops_zyg_poptable.csv"
vcf_dir="~/popgen/not_European_wine/"

##USAGE: modify the input paths, then run in the folder above sorted population folders

def parse_poptable(poptable_file):
    '''parse a csv file with
        -column of Standardized name as A
        -colum of Clades as B
        data from 1011 yeast genomes paper
        remove non-alphanumeric characters'''
    het_pop_dict={}
    with open(poptable_file) as f:
        next(f)
        for line in f:
            row_data=line.split(",") 
            strain=row_data[0]
            zygosity=row_data[1]
            pop=row_data[2].split("(")[0].strip()
            pop=pop.replace(" ","")
            pop=pop.replace(".","")
            pop=pop.replace("/","")

            if zygosity=="heterozygous":
                if pop=="":
                    pop="no_pop_assigned"
                if pop in het_pop_dict:
                    het_pop_dict[pop].append(strain)
                else:
                    het_pop_dict[pop]=[strain]
    return het_pop_dict

def sort_pops(het_pop_dict):
    #start in the vcf_dir file
    cmd="cd "+vcf_dir
    os.system(cmd)
    
    for pop in het_pop_dict:
        cmd = "cd "+vcf_dir+pop
        os.system(cmd)
        hetfolder=vcf_dir+pop+"/"+"heterozygotes"+"_"+pop
        
        if os.path.isdir(hetfolder)!=True:
            #make het  directory if it doesn't exist
            cmd = 'mkdir '+hetfolder
            os.system(cmd)
        for strain in het_pop_dict[pop]:
            cmd = 'cp '+vcf_dir+pop+"/"+strain+'.vcf* '+hetfolder+'/'
            #print(cmd)
            #exit()
            os.system(cmd)
            
            
het_pop_dict=parse_poptable(poptable_file)
sort_pops(het_pop_dict)



