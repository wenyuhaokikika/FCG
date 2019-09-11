# -*- coding:utf-8 -*-
#
# author: wenyuhao
# file: catch.py
# time: 2019/9/11
import pandas as pd
import numpy as np
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import os
import subprocess
import click
import re
from collections import defaultdict
import json
import multiprocessing as mp
import matplotlib.pyplot as plt

AUGUSTUS='augustus'#这里写augustus的安装路径。
d=defaultdict(bool)
#catch(0,'result.fa',gfffile,fnafile)
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
def run_genemark(fnafile):
	pass
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
def run_augustus(fnafile):
	Args=[AUGUSTUS,'--gff3=on','--outfile=augustus_out.gff3','--species=saccharomyces_cerevisiae_S288C',fnafile]
	p = subprocess.Popen(Args)#augustus --gff3=on --outfile=augustus_out1.gff3 --species=saccharomyces_cerevisiae_S288C out_little.fasta
	p.wait()#父进程等待子进程
	if p.poll() or p.poll()==0:#判断是否运行成功.
		print('augustus successfull exculate!!!')
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
def update_gff(gfffile,step):
	data=pd.DataFrame()
	if(step==0):
		data=pd.read_csv(gfffile,sep="\t",skiprows=1,header=None)
	else:
		data=pd.read_csv(gfffile,sep="\t",header=None)
	s=data[8]+";chrom="+data[0]+";start="+data[3].astype(str)+';end='+data[4].astype(str)+';strand='+data[6]
	data.loc[~s.apply(lambda x:d[x]),:].to_csv(gfffile,header=False,index=None,sep="\t")
def plotResult():
	a=[]
	with open('result/log.txt','r') as f:
		for line in f.readlines():
			if('>' in line):
				a.append(line)
	X,Y=[],[]
	for x,y in list(map(lambda x:x.replace('>','').replace(';number:',',').strip().split(','),a)):
		X.append(int(x))
		Y.append(int(y))
	fig, ax = plt.subplots(figsize=(18, 5))
	ax.plot(X, Y,'b')
	ax.set(xlabel='the length of the step', ylabel='the number of found complete gene',
			title='statics of iteration')
	ax.grid()
	fig.savefig("statics.png")
class FileNotFoundException(Exception):#这个异常就是文件找不到异常
	pass
class stepNotMatchNumsException(Exception):#这个异常指每次迭代的数目不是最终延长长度的倍数
	pass
@click.command()
@click.option('--step', default=1, help='每次迭代向两端衍生的长度。')
@click.option('--maxlen', default=10, help='向两端衍生的最长长度，必须是step的倍数，否侧抛出异常。')
@click.option('--gfffile', help='blast后perl的转化gff文件')
@click.option('--fnafile', help='原始基因组的fna文件')
@click.option('--Au',default='augustus',help='augustus的安装地址，默认是augustus，如果你加入到环境变量中，就不需要改变。')
def Iteration(step,gfffile,fnafile):
	'''
	这个程序切去衍生的序列做成gff文件，建模，\n
	得倒的文件筛选含有开始和结尾的文件，如果有开始和结尾，就认为这个预测基因完整。\n
	完整的基因记下后去除，不完整的进行下一轮迭代。\n
	/******contribute the code for projction if you are my brother,**********/\n
	/******The slacking is not my brother.***********************************/\n
	'''
	if(maxlen%step!=0):
		raise stepNotMatchNumsException('maxlen必须是step的倍数')
	global AUGUSTUS
	AUGUSTUS=Au
	run(gfffile,fnafile,step,maxlen)
def run(gfffile,fnafile,step,maxlen):
	if(os.path.exists('result')):
		os.system('rm -rf result')
	if not os.path.exists(gfffile):
		raise FileNotFoundException('blasted gff not found')
	os.system('mkdir result')
	os.system('cp {} result/'.format(gfffile))
	gfffile='result/'+gfffile
	resultfile='result.fa'#结果文件路径
	for step in range(0,maxlen,step):
		if(os.path.exists('result.fa')):
			os.system('rm -rf result.fa')
		catch(step,resultfile,gfffile,fnafile)
		if(os.path.exists('augustus_out.gff3')):
			os.system('rm -rf augustus_out.gff3')
		run_augustus(resultfile)
		pipline_augustus('augustus_out.gff3',step)
		update_gff(gfffile,step)
		with open('result/log.txt','a+') as f:
			f.write(">{};number:{}\n".format(str(step),str(len([k+'\n' for k,v in d.items() if v])))+"".join([k+'\n' for k,v in d.items() if v]))
	with open('result/finall.txt','a+') as f:
		f.write(json.dumps(d))
	plotResult()
if __name__=='__main__':
	Iteration()

