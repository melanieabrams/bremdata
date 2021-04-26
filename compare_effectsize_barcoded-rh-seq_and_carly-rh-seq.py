import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

#PARAMS

sep='\t'

unbarcoded_allfx='Carly_Reanalysis_1.5_1.1_coeffvarANDreadcutoff_uncollapsedMBA_unpaired_ratiocvfilt_20.0_2.0.all_fxsizes'
barcoded_37_allfx='37C_2.0_1.1_coeffvarANDreadcutoff_uncollapsedMBA_unpaired_ratiocvfilt_10.0_3.0.all_fxsizes'

goi=['YGR198W', 'YMR207C', 'YGL082W', 'YNL049C', 'YDL035C', 'YDR508C', 'YBR136W',
 'YML099C', 'YPL254W', 'YIL152W', 'YKL017C', 'YGR140W', 'YJR127C', 'YDR375C',
 'YOR091W', 'YLR397C', 'YNL132W', 'YMR078C', 'YLR422W', 'YMR125W', 'YOR371C',
 'YMR094W', 'YMR167W', 'YDR103W', 'YDR318W', 'YAL026C', 'YDR180W', 'YOR092W',
 'YDR235W', 'YER151C', 'YMR275C', 'YKL114C', 'YOL081W', 'YPR049C', 'YGL095C',
 'YDR456W', 'YKL197C', 'YIL068C', 'YOR326W', 'YNR045W', 'YJR107W', 'YPL268W',
 'YJL062W', 'YCR042C'] #hits 37_2.0_1.1_..._10.0_3.0 with sc defect and effect size >0.5


##
###goi: 2_1.1_1.6_3.0, pcutoff 0.05, abs(sp_mean_obs)<=0.75, sc_mean_obs<-0.5
##goi=['YGR198W','YLR397C','YBR136W','YGR140W','YNL049C','YKL017C','YGL082W',
##     'YMR125W','YDR508C','YMR207C','YDR375C','YDR180W','YKL197C','YDR318W',
##     'YMR094W','YOR326W','YBR081C','YPR049C','YIL152W','YER151C','YJR107W',
##     'YAL026C','YDR456W','YLR141W','YPL268W','YDR235W']

#RUN

ubcd=pd.read_csv(unbarcoded_allfx,sep=sep)
#print(ubcd.head())
ubcd['ubcd_nonabsolute_fx']=ubcd['sc_mean']-ubcd['sp_mean']
ubcd=ubcd.drop(columns=['effect_size','sc_mean','sp_mean','mwu_pval'])
ubcd=ubcd.set_index('gene')

bcd=pd.read_csv(barcoded_37_allfx,sep=sep)
bcd['bcd_nonabsolute_fx']=bcd['sc_mean']-bcd['sp_mean']
bcd=bcd.drop(columns=['effect_size','sc_mean','sp_mean','mwu_pval'])
bcd=bcd.set_index('gene')

both=pd.concat([ubcd,bcd],axis=1, join="inner")
#print(both.head())

both=both.reset_index()
#print(both.head())
both_goi=both[both['gene'].isin(goi)]
print(both_goi)

fig, ax = plt.subplots(figsize=(6,6))
plt.scatter(both['bcd_nonabsolute_fx'],both['ubcd_nonabsolute_fx'],color='black',s=5)
plt.scatter(both_goi['bcd_nonabsolute_fx'],both_goi['ubcd_nonabsolute_fx'],color='red',s=10)
plt.xlabel('37C_barcoded fxsize')
plt.ylabel('39C unbarcoded fxsize')
plt.ylim(-3,3)
plt.xlim(-3,3)
plt.hlines(0,-3,3,colors='black')
plt.show()



