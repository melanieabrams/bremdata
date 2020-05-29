import sys

##USAGE: python alignments_get_only_spar.py alignment1.fa alignment2.fa
#assumes alignments named like Scer_Spar_YGR098C_YGR098C.ga

popname="1WineEuropean"
strains=['BFC','AGL','ABE','SACE_YAM','BKI','ABR','ABQ','BNV','BPB','BPT','BQM','BFK','BQS','BQE','AIN','BQR','BLT','AIF','AQR','BPF','BLR','BNL','BPK','BKQ','BPG','BPH','BQT','APT','BID','AGI','BHV','BNE','BLA','BKT','BIA','BKL','BIB','BKK','BLQ','SACE_YBY','SACE_YDF','SACE_YDO','AHG','AII','AHD','AGF','AHE','AQQ','ASG','ASD','ASE','AIQ','BLG','BLI','ACP','BQA','AMA','ALQ','BBI','ALT','BBK','BFB','BFD','ALS','BBH','BPQ','BQN','CBI','CCM','BEP','BRC','CGT','ALL','ALM','AVH','AHR','SACE_YAD','AGH','BTH','ATH','BFQ','CIA','BFR','CIB','CBG','ACF','SACE_YCD','APS','ARC','BER','BES','ARK','ATA','ARQ','SACE_YBN','APV','CFV','SACE_YDC','SACE_GAP','AAI','ARV','BHC','BHF','BLB','BPP','BTN','BTL','CRH','SACE_YAA','BSP','AHT','ATN','SACE_YBO','BNK','BNR','BQH','BEQ','AIK','CAH','CBA','CAI','CAK','CAL','BSL','SACE_YBX','BSE','AAA','CFS','BEV','BET','BST','ASS','BSD','ATM','CQT','BSA','ADR','BSC','CFQ','CFR','BSB','ADC','SACE_YCP','SACE_YCQ','BII','BIS','BIG','AIR','AIT','BNF','BAR','BAS','BAT','BRD','AGV','BQI','AGA','AGS','BPC','BKE','BLH','CRG','AIL','AMS','BLE','ADB','BQB','AMR','ABT','BRA','BIC','BLK','BBV','BVK','CFE','BRF','BSM','BSN','CPE','BLM','BQC','BRE','BPV','CRC','BHN','BBB','BAV','CAE','CAF','CRL','BBA','SACE_YAY','ADI','SACE_YBA','SACE_YAZ','BMI','SACE_YBE','BMK','SACE_YBF','BMH','ADD','SACE_YBB','SACE_YBG','SACE_YBH','CDS','BIM','BIN','BKA','BQQ','ADS','BNM','ACR','BIH','BKC','BLS','BKD','BNP','BKN','BNQ','CBH','BKH','BPM','BKP','CRA','CRB','CHC','CGS','CHB','CHD','CGV','CHA','CHE','AHI','BKM','BIF','BKG','BKF','BLL','ARL','CKI','CKL','CKH','CKD','CKE','CIS','BMN','CIT','CKP','CKQ','CKB','CKC','CIQ','CIM','CIN','CIR','CKF','CKG','CKM','CIP','BCA','AMT','CIV','CKT','BIL','ABS','BIP','BIQ','BNI','BLV','BQG','CRD','CRI','BVM','CQV','CRK','BLF','CEA','CEB','BNT','BLD','BKR','CEC','BNH','BLC','CEE','BNN','CED','BKS','BLP','BLN','CDV','CEF','ANQ','AVK','BIK','CLH','AHN','BSG','SACE_YCA','ACT','ACV','AKP','AKV','ALE','BQP','BRB','BIT','AQB','BRH','ABP','BIR','CPF','BNG','AIM','CFB','CHF','BHB','AAK','AHP','AET','AEK','AQE','AQP','ASA','BPL','AIP','CAN','CAB','CAP','AQC','AHK','AES','AFI','BPR','AIA','BVT','CAA','CAC','BVV','CAD','CAM','CAQ','CBB','CBE','AAE','AIB','CAS','CAR','CAT','CBC','CAV','CBD','AFV','ART','AIC','SACE_YAF','BSR','BSS','BTR','BTP','ASR','SACE_YDJ']

def get_only_onepop(alignment):
    savename=popname+(alignment.replace("Spar_","")).split('.')[0]+'.fasta'
    with open(savename, 'w') as wf:
        f=open(alignment)
        for line in f:
            if line[0]==">":
                writeseq=False
                if line[7:10] in strains: 
                    writeseq=True #write the sequence following this header if this is a scer strain in the list
                    wf.writelines('>'+line[7:10]+'\n')
            else:
                if writeseq==True:
                    wf.writelines(line)
                    
        f.close()

#RUN
alignments=sys.argv[1:]
for alignment in alignments:
    get_only_onepop(alignment)
