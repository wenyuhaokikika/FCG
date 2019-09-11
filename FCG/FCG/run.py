# -*- coding:utf-8 -*-
#
# author: wenyuhao
# file: run.py
# time: 2019/9/11

def run_augustus(fnafile):
	Args=[AUGUSTUS,'--gff3=on','--outfile=augustus_out.gff3','--species=saccharomyces_cerevisiae_S288C',fnafile]
	p = subprocess.Popen(Args)#augustus --gff3=on --outfile=augustus_out1.gff3 --species=saccharomyces_cerevisiae_S288C out_little.fasta
	p.wait()#父进程等待子进程
	if p.poll() or p.poll()==0:#判断是否运行成功.
		print('augustus successfull exculate!!!')
def run_genemark(fnafile):
	pass
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