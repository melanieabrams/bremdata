import sys

###USAGE

#python assignSampleNames_h12haps.py #all else in parameters

###PARAMETERS

#WINE EUROPEAN
#sample_order=['AAA', 'AAE', 'AAI', 'AAK', 'ABE', 'ABP', 'ABQ', 'ABR', 'ABS', 'ABT', 'ACF', 'ACP', 'ACR', 'ACT', 'ACV', 'ADB', 'ADC', 'ADD', 'ADI', 'ADR', 'ADS', 'AEK', 'AES', 'AET', 'AFI', 'AFV', 'AGA', 'AGF', 'AGH', 'AGI', 'AGL', 'AGS', 'AGV', 'AHD', 'AHE', 'AHG', 'AHI', 'AHK', 'AHN', 'AHP', 'AHR', 'AHT', 'AIA', 'AIB', 'AIC', 'AIF', 'AII', 'AIK', 'AIL', 'AIM', 'AIN', 'AIP', 'AIQ', 'AIR', 'AIT', 'AKP', 'AKV', 'ALE', 'ALL', 'ALM', 'ALQ', 'ALS', 'ALT', 'AMA', 'AMR', 'AMS', 'AMT', 'ANQ', 'APS', 'APT', 'APV', 'AQB', 'AQC', 'AQE', 'AQP', 'AQQ', 'AQR', 'ARC', 'ARK', 'ARL', 'ARQ', 'ART', 'ARV', 'ASA', 'ASD', 'ASE', 'ASG', 'ASR', 'ASS', 'ATA', 'ATH', 'ATM', 'ATN', 'AVH', 'AVK', 'BAR', 'BAS', 'BAT', 'BAV', 'BBA', 'BBB', 'BBH', 'BBI', 'BBK', 'BBV', 'BCA', 'BEP', 'BEQ', 'BER', 'BES', 'BET', 'BEV', 'BFB', 'BFC', 'BFD', 'BFK', 'BFQ', 'BFR', 'BHB', 'BHC', 'BHF', 'BHN', 'BHV', 'BIA', 'BIB', 'BIC', 'BID', 'BIF', 'BIG', 'BIH', 'BII', 'BIK', 'BIL', 'BIM', 'BIN', 'BIP', 'BIQ', 'BIR', 'BIS', 'BIT', 'BKA', 'BKC', 'BKD', 'BKE', 'BKF', 'BKG', 'BKH', 'BKI', 'BKK', 'BKL', 'BKM', 'BKN', 'BKP', 'BKQ', 'BKR', 'BKS', 'BKT', 'BLA', 'BLB', 'BLC', 'BLD', 'BLE', 'BLF', 'BLG', 'BLH', 'BLI', 'BLK', 'BLL', 'BLM', 'BLN', 'BLP', 'BLQ', 'BLR', 'BLS', 'BLT', 'BLV', 'BMH', 'BMI', 'BMK', 'BMN', 'BNE', 'BNF', 'BNG', 'BNH', 'BNI', 'BNK', 'BNL', 'BNM', 'BNN', 'BNP', 'BNQ', 'BNR', 'BNT', 'BNV', 'BPB', 'BPC', 'BPF', 'BPG', 'BPH', 'BPK', 'BPL', 'BPM', 'BPP', 'BPQ', 'BPR', 'BPT', 'BPV', 'BQA', 'BQB', 'BQC', 'BQE', 'BQG', 'BQH', 'BQI', 'BQM', 'BQN', 'BQP', 'BQQ', 'BQR', 'BQS', 'BQT', 'BRA', 'BRB', 'BRC', 'BRD', 'BRE', 'BRF', 'BRH', 'BSA', 'BSB', 'BSC', 'BSD', 'BSE', 'BSG', 'BSL', 'BSM', 'BSN', 'BSP', 'BSR', 'BSS', 'BST', 'BTH', 'BTL', 'BTN', 'BTP', 'BTR', 'BVK', 'BVM', 'BVT', 'BVV', 'CAA', 'CAB', 'CAC', 'CAD', 'CAE', 'CAF', 'CAH', 'CAI', 'CAK', 'CAL', 'CAM', 'CAN', 'CAP', 'CAQ', 'CAR', 'CAS', 'CAT', 'CAV', 'CBA', 'CBB', 'CBC', 'CBD', 'CBE', 'CBG', 'CBH', 'CBI', 'CCM', 'CDS', 'CDV', 'CEA', 'CEB', 'CEC', 'CED', 'CEE', 'CEF', 'CFB', 'CFE', 'CFQ', 'CFR', 'CFS', 'CFV', 'CGS', 'CGT', 'CGV', 'CHA', 'CHB', 'CHC', 'CHD', 'CHE', 'CHF', 'CIA', 'CIB', 'CIM', 'CIN', 'CIP', 'CIQ', 'CIR', 'CIS', 'CIT', 'CIV', 'CKB', 'CKC', 'CKD', 'CKE', 'CKF', 'CKG', 'CKH', 'CKI', 'CKL', 'CKM', 'CKP', 'CKQ', 'CKT', 'CLH', 'CPE', 'CPF', 'CQT', 'CQV', 'CRA', 'CRB', 'CRC', 'CRD', 'CRG', 'CRH', 'CRI', 'CRK', 'CRL', 'SACE_GAP', 'SACE_YAA', 'SACE_YAD', 'SACE_YAF', 'SACE_YAM', 'SACE_YAY', 'SACE_YAZ', 'SACE_YBA', 'SACE_YBB', 'SACE_YBE', 'SACE_YBF', 'SACE_YBG', 'SACE_YBH', 'SACE_YBN', 'SACE_YBO', 'SACE_YBX', 'SACE_YBY', 'SACE_YCA', 'SACE_YCD', 'SACE_YCP', 'SACE_YCQ', 'SACE_YDC', 'SACE_YDF', 'SACE_YDJ', 'SACE_YDO']
#sample order generated sample_order using vcf2SelectHapStats_getSampleOrderKey.py on head (1000) of vcf which replaces the allele or '.' with the sample name in the vcf2shsinput generation
#h12_h2h1='merged_1WineEuropean_chromosome7.h12_h2h1'

#BRAZILIAN BIOETHANOL
sample_order =['AEG', 'AEH', 'AFD', 'AFR', 'AGM', 'AHC', 'APG', 'APH', 'BTT', 'BTV', 'BVA', 'BVB', 'BVC', 'BVD', 'BVE', 'BVF', 'BVG', 'BVH', 'CMT', 'CNB', 'CNE', 'CNF', 'CNG', 'CNH', 'CNI', 'CNK', 'CNL', 'CNM', 'CNN', 'CNP', 'CNQ', 'CNR', 'CNS', 'CNT', 'CNV']
h12_h2h1='merged_3BrazilianBioethanol_chromosome7.h12_h2h1'

#GOI
gene='YGR098C'

#GFF
gff='/usr2/people/mabrams/Amended_Genomes/S288C/saccharomyces_cerevisiae_R64-1-1_20110208_withoutFasta.gff'
roman_numerals_in_gff=True




###

def ParseFromGFF(gfffile):
    '''
    Parses SGD features flat file
    Input: SGD_features.tab file
    Output: dict of {chrom:{gene:[start,stop]}}
    '''
    roman_to_numerals={
        'chrI':'chr01','chrII':'chr02','chrIII':'chr03','chrIV':'chr04','chrV':'chr05',
        'chrVI':'chr06','chrVII':'chr07','chrVIII':'chr08','chrIX':'chr09','chrX':'chr10',
        'chrXI':'chr11','chrXII':'chr12','chrXIII':'chr13','chrXIV':'chr14','chrXV':'chr15',
        'chrXVI':'chr16','chrMito':'chrMito','2-micron':'chr2u'}
    gff_genes={}
    
    f = open(gfffile)
    lines=[]
    for line in f:
        if line[0]!='#' : #skip header rows
            row_data=line.split("\t")
            chrom=row_data[0]
            start=int(row_data[3])
            stop=int(row_data[4])
            gene_length=stop - start + 1
            info=row_data[8].split(";")
            yName=info[0].split('=')[1]
            if yName[0]=='Y' and len(yName)>5:
                gff_genes[yName]=(start,stop)
                
    f.close()
    
    return gff_genes



def AmendH12(h12_file):
    '''
    replaces strain number with strain name into a new outfile
    '''

    with open(outfile_name,'w') as wf:
        #add position and values for each base to a dictionary for that chromosome
        f = open(h12_file)
        next(f)
        for line in f:
            row_data = line.strip().split("\t")
            center=int(row_data[0])
            if gff_genes[gene][0]<center<gff_genes[gene][1]: #only pull points where center window was in gene

                strainnums=row_data[5]
                split_strainnums=strainnums.replace('(','').split(')')
                replaced=''
                for i in range(len(split_strainnums)):
                    hapstrains=''
                    hapstrainnums=split_strainnums[i].split(',')
                    for strainnum in hapstrainnums:
                        if strainnum!='':
                            hapstrains+=sample_order[int(strainnum)-1]+','
                    #print(hapstrains)
                    replaced+='('+hapstrains[:-1]+')'
                newline=''
                for i in range(len(row_data)):
                    if i!=5:
                        newline+=row_data[i]+'\t'
                    else:
                        newline+=replaced+'\t'
                wf.writelines(newline+'\n')
            
        
        #columns:
            #1ctrcoord (index 0)
            #2leftcoord
            #3rightcoord
            #4K
            #55hapfreqspec
            #6strainnum
            #7H1 (index 6)
            #8H2 (index 7)
            #9H12 (index 8)
            #10H2/H1 (index 9)
            #11H123 (index 10)
 
    f.close()



    
    



##RUN###
outfile_name=gene+'_'+h12_h2h1.split('.')[0]+'.h12_h2h1'
gff_genes=ParseFromGFF(gff)
AmendH12(h12_h2h1)
