# -*- coding:utf-8 -*-
#
# author: wenyuhao
# file: update.py
# time: 2019/9/11
import pandas as pd
def update_gff(gfffile,step):
	data=pd.DataFrame()
	if(step==0):
		data=pd.read_csv(gfffile,sep="\t",skiprows=1,header=None)
	else:
		data=pd.read_csv(gfffile,sep="\t",header=None)
	s=data[8]+";chrom="+data[0]+";start="+data[3].astype(str)+';end='+data[4].astype(str)+';strand='+data[6]
	data.loc[~s.apply(lambda x:d[x]),:].to_csv(gfffile,header=False,index=None,sep="\t")