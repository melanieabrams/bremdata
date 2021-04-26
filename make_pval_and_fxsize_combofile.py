import pandas as pd

#PARAMS
sep='\t'
prefix='37C_2.0_1.1_coeffvarANDreadcutoff_uncollapsedMBA_unpaired_ratiocvfilt_1.6_3.0'

#RUN
mwu_results=prefix+'.mwu_test_results'
fxsize_results=prefix+'.fxsizes'

mwu_df=pd.read_csv(mwu_results,sep=sep)
fxsize_df = pd.read_csv(fxsize_results,sep=sep)


both=fxsize_df.merge(mwu_df)
both=both.set_index('gene')
both.to_csv(prefix+'.combined_p_and_fx',sep=sep)
