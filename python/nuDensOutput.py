"""Functions and classes that read the output and help to do plots"""

import re
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import interp1d

colors = ["r", "g", "b", "k", "c", "m", "y", "#cc6699", "#ff9933"]
styles = ["-", "--", ":", "-."]

def finalizePlot(fname,
		lloc="best",
		xlab=None,
		ylab=None,
		xscale=None,
		yscale=None,
		xlim=None,
		ylim=None,
		):
	ax=plt.gca()
	ax.tick_params("both", which="both", direction="out",
		left=True, right=True, top=True, bottom=True)
	plt.legend(loc=lloc)
	if xlab is not None:
		plt.xlabel(xlab)
	if ylab is not None:
		plt.ylabel(ylab)
	if xscale is not None:
		plt.xscale(xscale)
	if yscale is not None:
		plt.yscale(yscale)
	if xlim is not None:
		plt.xlim(xlim)
	if ylim is not None:
		plt.ylim(ylim)
	plt.tight_layout()
	plt.savefig(fname)
	plt.close()

def stripRepeated(data, ix1, ix2):
	"""Strip the repeated points from an output file,
	to avoid steps in the plots
	"""
	x=[]
	y=[]
	for d in data:
		if len(y) == 0:
			x.append(d[ix1])
			y.append(d[ix2])
		if y[-1] != d[ix2]:
			x.append(d[ix1])
			y.append(d[ix2])
	x.append(data[-1][ix1])
	y.append(data[-1][ix2])
	return np.asarray(x), np.asarray(y)

class NuDensRun():
	"""Class that reads the output and helps to do plots"""

	def __init__(self, folder, nnu=3, full=True, label=""):
		self.folder = folder
		self.full = full
		self.label = label
		fdy = np.loadtxt("%s/fd.dat"%folder)
		self.yv = fdy[:,0]
		self.fd = fdy[:,1]
		self.zdat = np.loadtxt("%s/z.dat"%folder)
		self.nnu = nnu
		self.rho = np.asarray([
			[[None, None] for i in range(nnu)] for j in range(nnu)
			])
		try:
			with open("%s/resume.dat"%folder) as _f:
				self.resume = _f.readlines()
		except FileNotFoundError:
			self.resume = [""]*(self.nnu+2)
			self.hasResume= False
		else:
			self.hasResume = True
		if self.hasResume:
			self.Neff = float(
				re.match("Neff[ =]*([\d.]*)", self.resume[-1]).group(1))
			self.zfin = float(
				re.match("final z[ =]*([\d.]*)", self.resume[0]).group(1))
		self.deltarhofin = []
		for i in range(self.nnu):
			if self.hasResume:
				self.deltarhofin.append(
					float(
						re.match("dRho_%s[ =]*([\d.]*)"%(i+1), self.resume[i+1])
							.group(1)
						)
					)
			self.rho[i, i, 0] = np.loadtxt(
				"%s/nuDens_diag%d.dat"%(folder, i+1))
			if full:
				for j in range(i+1, self.nnu):
					self.rho[i, j, 0] = np.loadtxt(
						"%s/nuDens_nd_%d%d_re.dat"%(folder, i+1, j+1))
					self.rho[i, j, 1] = np.loadtxt(
						"%s/nuDens_nd_%d%d_im.dat"%(folder, i+1, j+1))
		self.printTableLine()

	def printTableLine(self):
		if self.hasResume:
			deltastr = ""
			for i in range(self.nnu):
				deltastr += "{:.5f} & ".format(self.deltarhofin[i])
			print("{lab:<15s} & {zfin:.5f} & {deltastr:s}{Neff:.5f}\\\\".format(
				lab=self.label,
				Neff=self.Neff,
				zfin=self.zfin,
				deltastr=deltastr,
				))

	def plotFD(self, ls="-", lc="k"):
		plt.plot(
			self.yv, self.fd,
			label=self.label, ls=ls, marker=".", c=lc
			)
		plt.xscale("log")
		plt.yscale("log")
		plt.xlabel("$y$")
		plt.ylabel(r"$y^2 f(y)$")

	def plotZ(self, ls="-", lc="k"):
		plt.plot(
			*stripRepeated(self.zdat, 0, 1),
			label=self.label, ls=ls, c=lc
			)
		plt.xscale("log")
		plt.xlabel("$x$")
		plt.ylabel(r"$z$")

	def plotDeltaZ(self, ref, ls="-", lc="k"):
		mex, mey = stripRepeated(self.zdat, 0, 1)
		mef = interp1d(mex, mey)
		refx, refy = stripRepeated(ref.zdat, 0, 1)
		reff = interp1d(refx, refy)
		plt.plot(
			mex, reff(mex)-mef(mex),
			label=self.label, ls=ls, c=lc
			)
		plt.xscale("log")
		plt.xlabel("$x$")
		plt.ylabel(r"$z$")

	def plotRhoDiag(self, inu, iy, ls, lc="k"):
		plt.plot(
			*stripRepeated(self.rho[inu, inu, 0], 0, iy),
			label="%s %d"%(self.label, inu+1), ls=ls, c=lc
			)
		plt.xscale("log")
		plt.xlabel("$x$")
		plt.ylabel(r"$\rho/\rho_{eq}-1$")

	def plotRhoOffDiag(self, i1, i2, iy, lc="k"):
		if not self.full:
			print("no offdiagonal loaded")
			return
		plt.plot(
			*stripRepeated(self.rho[i1, i2, 0], 0, iy),
			ls="-", c=lc, label="%s %d%d re"%(self.label, i1+1, i2+1)
			)
		plt.plot(
			*stripRepeated(self.rho[i1, i2, 1], 0, iy),
			ls=":", c=lc, label="%s %d%d im"%(self.label, i1+1, i2+1)
			)
		plt.xscale("log")
		plt.xlabel("$x$")
		plt.ylabel(r"$\rho_{ij}$")

	def plotdRhoOffDiag(self, i1, i2, iy, lc="k"):
		if not self.full:
			print("no offdiagonal loaded")
			return
		dijrex, dijrey = stripRepeated(self.rho[i1, i2, 0], 0, iy)
		dijimx, dijimy = stripRepeated(self.rho[i1, i2, 1], 0, iy)
		plt.plot(
			dijrex, np.gradient(dijrey, dijrex),
			ls="-", c=lc, label="%s %d%d re"%(self.label, i1+1, i2+1)
			)
		plt.plot(
			dijimx, np.gradient(dijimy, dijimx),
			ls=":", c=lc, label="%s %d%d im"%(self.label, i1+1, i2+1)
			)
		plt.xscale("log")
		plt.xlabel("$x$")
		plt.ylabel(r"$\rho_{ij}$")

	def plotRhoFin(self, ix, iy=None, ri=0, ls="-", lc="k"):
		if iy is None:
			iy = ix
		if ri not in [0, 1]:
			ri = 0
		plt.plot(
			self.yv, self.rho[ix, iy, ri][-1, 1:],
			ls=ls, c=lc,
			label="%s %d%d %s"%(self.label, ix+1, iy+1, "re" if ri == 0 else "im"),
			)
		plt.xlabel("$y$")
		plt.ylabel(r"$\rho_{ij}^{\rm fin}(y)$")
