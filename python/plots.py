import matplotlib, sys, time
#matplotlib.use('Agg')
import numpy as np
import matplotlib.pyplot as plt

from matplotlib.backends.backend_pdf import PdfPages
pdf = PdfPages('plots/last.pdf')
ymin=0.01
ymax=20
ny=100
yarr=np.logspace(np.log10(ymin),np.log10(ymax),ny)

columns=[0,50,-1]
def plot_momentum(fname, columns):
	pts=np.loadtxt(fname)
	ptsStruc=[]
	for line in pts:
		tmp=line[1:] 
		ptsStruc.append([line[0], tmp])
		
	plt.subplot(111)
	plt.xscale('log')
	#plt.xlabel("weight")
	#plt.ylabel("loglikelihood")
	#plt.xlim([1e-10,1])
	#plt.ylim([11260,11450])
	for q in columns:
		print q, ptsStruc[q][0]#, ptsStruc[q][1]
		plt.plot(yarr,ptsStruc[q][1],label="%f"%ptsStruc[q][0])
	plt.legend(loc='upper right')
	pdf.savefig()#time.strftime("%y%m%d_%H%M%S")+
	#if not noShow:
	plt.show()
	plt.close()
def plot_z(fname):
	pts=np.loadtxt(fname)
	x=[]
	y=[]
	for line in pts:
		x.append(line[0]) 
		y.append(line[1])
		
	plt.subplot(111)
	plt.xscale('log')
	#plt.xlabel("weight")
	#plt.ylabel("loglikelihood")
	#plt.xlim([1e-10,1])
	#plt.ylim([11260,11450])
	plt.plot(x,y,label="z")
	#plt.legend(loc='upper right')
	pdf.savefig()#time.strftime("%y%m%d_%H%M%S")+
	#if not noShow:
	plt.show()
	plt.close()

plot_momentum("output/nuDens_diag1.dat", columns)
plot_momentum("output/nuDens_diag2.dat", columns)
plot_momentum("output/nuDens_diag3.dat", columns)
plot_z("output/z.dat")

pdf.close()
