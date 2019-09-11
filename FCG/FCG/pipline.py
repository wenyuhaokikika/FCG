# -*- coding:utf-8 -*-
#
# author: wenyuhao
# file: pipline.py
# time: 2019/9/11
import pandas as pd
def pipline_augustus(gfffile,step):
	s,start_codon,end_codon,g='',False,False,''
	def parseline(line):
		return re.search(r'[(](.*?)[)]',line).group(1).split(', ')[1].replace('name =','').strip()
	with open(gfffile,'r') as f:
		for line in f.readlines():
			if('prediction on sequence number' in line):
				g=parseline(line)
				if not d[g]:
					start_codon,end_codon=False,False
			if('start_codon' in line):
				start_codon=True
			if('stop_codon' in line):
				end_codon=True
			if(start_codon&end_codon):
				d[g]=start_codon&end_codon
	with open(gfffile,'r') as f:
		with open('result/{}.gff'.format(str(step)),'a+') as sf:
			pan=True#写入指针，有起始和终止的基因不写入
			for line in f.readlines():
				if(pan):
					sf.write(line)
				if('prediction on sequence number' in line and not d[parseline(line)]):
					pan=False
				if('prediction on sequence number' in line and d[parseline(line)]):
					pan=True
#pipline_augustus(gfffile,0)
def pipline_genemark(gfffile,step):
	def parseFirstColumns(line):
		#s=re.search(r'[(](.*?)[)]',line).group(1).split(', ')[1].replace('name =','').strip()
		change=map(lambda x:{x.split('=')[0]:x.split('=')[1]},line.split(';'))
		di={}
		for i in change:
			di.update(i)
		return di
	def antiParseFirstColumns(di):
		return "".join(['{}={};'.format(str(k),str(v)) for k,v in di.items()])
	data=pd.read_csv(gfffile,sep="\t",header=None)
	for chrom,da in data.groupby(data[0]):
		if('start_codon' in da[2] and 'stop_codon' in da[2]):
			d[chrom]=True
	data.loc[~data.apply(lambda x:d[x[0]],axis=1),:].to_csv(gfffile,sep="\t",index=False,header=None)
	data.loc[data.apply(lambda x:d[x[0]],axis=1),:].to_csv('result/{}.gff'.format(str(step)),sep="\t",index=False,header=None)
