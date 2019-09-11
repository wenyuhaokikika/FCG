# -*- coding:utf-8 -*-
#
# author: wenyuhao
# file: catch.py
# time: 2019/9/11
import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import os
def catch(step,resultfile,gfffile,fnafile):
	if(os.path.exists(resultfile)):
		os.remove(resultfile)
	data=pd.DataFrame()
	if(step==0):
		data=pd.read_csv(gfffile,sep="\t",skiprows=1,header=None)
	else:
		data=pd.read_csv(gfffile,sep="\t",header=None)
	data=data.loc[data[2]=='gene',:]
	data[8]=data[8]+";chrom="+data[0]+";start="+data[3].astype(str)+';end='+data[4].astype(str)+';strand='+data[6]
	with open(resultfile,'a+') as sf:
		for seq_record in SeqIO.parse(fnafile, "fasta"):
			da=data.loc[data[0]==seq_record.id,:]
			da[3],da[4]=da[3]-step,da[4]+step
			for start,end,des in zip(da[3],da[4],da[8]):
				seq=seq_record.seq[start:end+1]
				re=SeqRecord(Seq(str(seq)),id=des,description='')
				SeqIO.write(re, sf, "fasta")