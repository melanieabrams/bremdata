from Bio import SeqIO

###GNP1/YDR508C
##infile='OG1243_YDR508C.mfa'
##outfile='YDR508C_pos255and256_context.mfa'
##first_codon=255
##last_codon=256
##
###CBF2/YGR140W
##infile='OG1967_YGR140W.mfa'
##outfile='YGR140W_pos505_context.mfa'
##first_codon=533 #because of cleandata parameter in PAML vs alignment
##last_codon=534


#CTF18/YMR078C
infile='OG3779_YMR078C'
outfile='YMR078C_pos286_context.mfa'
first_codon=286
last_codon=287




#PARAMETERS
rename_dict={'>Smik_':'>Smikatae',
             '>Spar_':'>Sparadoxus',
             '>Scer_':'>Scerevisiae',
             '>Skud_':'>Skudriavzevii',
             '>Sbay_':'>Suvarum'}
    
context_aa=30 #number of aa to include on each side

##FUNCTIONS###

def getRename(headerline):
    for item in rename_dict:
        if headerline.startswith(item):
            return rename_dict[item]
    
def parseAlignment(infile):
    alignment_dict={}
    sequence=''
    with open(infile) as f:
        for line in f:
            if line[0]=='>':
                #add sequence to alignment dict if there is one
                if sequence!='':
                    alignment_dict[seqname]=sequence
                #start new item 
                seqname=getRename(line) # get intelligable name
                sequence='' #reset sequence
            else:
                sequence+=line.strip()
        alignment_dict[seqname]=sequence #add last one
    return alignment_dict

def writeContext(alignment_dict,first_codon,last_codon,context_aa,outfile):
    start_index=3*(first_codon-context_aa)
    stop_index=3*(last_codon+context_aa)

    print(start_index)
    print(stop_index)

    
    with open(outfile,'w') as wf:
        for seqname in alignment_dict:
            print(seqname)
            wf.writelines(seqname+'\n')
            wf.writelines(alignment_dict[seqname][start_index:stop_index]+'\n')
    return
            
    

##RUN###
alignment_dict=parseAlignment(infile)

writeContext(alignment_dict,first_codon,last_codon,context_aa,outfile)
