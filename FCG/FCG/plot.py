# -*- coding:utf-8 -*-
#
# author: wenyuhao
# file: plot.py
# time: 2019/9/11
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