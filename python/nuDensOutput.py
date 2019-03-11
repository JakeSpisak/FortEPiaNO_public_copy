"""Functions and classes that read the output and help to do plots"""

import re
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import interp1d

colors = ["r", "g", "b", "k", "c", "m", "y", "#cc6699", "#ff9933"]
styles = ["-", "--", ":", "-."]

def finalizePlot(fname,
		lloc="best",
		title="",
		xlab=None,
		ylab=None,
		xscale=None,
		yscale=None,
		xlim=None,
		ylim=None,
		):
	plt.title(title)
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

	def __init__(self, folder, nnu=3, full=True, label="", plots=False):
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
				re.match("Neff[ =]*([-\d.]*)", self.resume[-1]).group(1))
			self.zfin = float(
				re.match("final z[ =]*([-\d.]*)", self.resume[0]).group(1))
		self.deltarhofin = []
		for i in range(self.nnu):
			if self.hasResume:
				self.deltarhofin.append(
					float(
						re.match("dRho_%s[ =]*([-\d.]*)"%(i+1), self.resume[i+1])
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
		if plots:
			self.doAllPlots()

	def interpolateRhoIJ(self, i1, i2, y, ri=0):
		xv = []
		yv = []
		prevy = 0
		for i, x in enumerate(self.rho[i1, i2, ri][:, 0]):
			fy = interp1d(self.yv, self.rho[i1, i2, ri][i, 1:])
			cy = fy(y)
			if cy != prevy:
				prevy = cy
				yv.append(prevy)
				xv.append(x)
		xv.append(x)
		yv.append(cy)
		return xv, yv

	def interpolateRhoIJ_x(self, i1, i2, x, ri=0):
		ov = []
		for i, y in enumerate(self.yv):
			fx = interp1d(self.rho[i1, i2, ri][:, 0], self.rho[i1, i2, ri][:, i])
			ov.append(fx(x))
		return self.yv, ov

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
		plt.ylabel(r"$z-z_{\rm ref}$")

	def plotRhoDiag(self, inu, iy, ls, lc="k"):
		plt.plot(
			*stripRepeated(self.rho[inu, inu, 0], 0, iy),
			label="%s %d"%(self.label, inu+1), ls=ls, c=lc
			)
		plt.xscale("log")
		plt.xlabel("$x$")
		plt.ylabel(r"$\rho/\rho_{eq}-1$")

	def plotRhoOffDiag(self, i1, i2, iy, lc="k", im=True):
		if not self.full:
			print("no offdiagonal loaded")
			return
		plt.plot(
			*stripRepeated(self.rho[i1, i2, 0], 0, iy),
			ls="-", c=lc, label="%s %d%d re"%(self.label, i1+1, i2+1)
			)
		if im:
			plt.plot(
				*stripRepeated(self.rho[i1, i2, 1], 0, iy),
				ls=":", c=lc, label="%s %d%d im"%(self.label, i1+1, i2+1)
				)
		plt.xscale("log")
		plt.xlabel("$x$")
		plt.ylabel(r"$\rho_{ij}$")

	def plotdRhoOffDiag(self, i1, i2, iy, lc="k", im=True):
		if not self.full:
			print("no offdiagonal loaded")
			return
		dijrex, dijrey = stripRepeated(self.rho[i1, i2, 0], 0, iy)
		plt.plot(
			dijrex, np.gradient(dijrey, dijrex),
			ls="-", c=lc, label="%s %d%d re"%(self.label, i1+1, i2+1)
			)
		if im:
			dijimx, dijimy = stripRepeated(self.rho[i1, i2, 1], 0, iy)
			plt.plot(
				dijimx, np.gradient(dijimy, dijimx),
				ls=":", c=lc, label="%s %d%d im"%(self.label, i1+1, i2+1)
				)
		plt.xscale("log")
		plt.xlabel("$x$")
		plt.ylabel(r"$d\rho_{ij}/dt$")

	def plotRhoFin(self, ix, iy=None, ri=0, ls="-", lc="k"):
		ylabel = r"$\rho_{ij}^{\rm fin}(y)$"
		if iy is None:
			iy = ix
			ylabel = r"$\rho_{ij}^{\rm fin}(y)/\rho_{eq}-1$"
		if ri not in [0, 1]:
			ri = 0
		plt.plot(
			self.yv, self.rho[ix, iy, ri][-1, 1:],
			ls=ls, c=lc,
			label="%s %d%d %s"%(self.label, ix+1, iy+1, "re" if ri == 0 else "im"),
			)
		plt.xlabel("$y$")
		plt.ylabel(ylabel)

	def plotRhoX(self, i1, x, i2=None, ri=0, ls="-", lc="k"):
		ylabel = r"$\rho_{ij}(y)$"
		if i2 is None:
			i2 = i1
			ylabel = r"$\rho_{ij}(y)/\rho_{eq}-1$"
		if ri not in [0, 1]:
			ri = 0
		plt.plot(
			*self.interpolateRhoIJ_x(i1, i2, x, ri),
			ls=ls, c=lc,
			label="%s %d%d %s x=%f"%(self.label, i1+1, i2+1, "re" if ri == 0 else "im", x),
			)
		plt.xlabel("$y$")
		plt.ylabel(ylabel)

	def plotRhoDiagY(self, inu, y, ls, lc="k"):
		plt.plot(
			*self.interpolateRhoIJ(inu, inu, y, ri=0),
			label="%s %d"%(self.label, inu+1), ls=ls, c=lc
			)
		plt.xscale("log")
		plt.xlabel("$x$")
		plt.ylabel(r"$\rho/\rho_{eq}-1$")

	def plotRhoOffDiagY(self, i1, i2, y, lc="k", ls="-", im=True):
		if not self.full:
			print("no offdiagonal loaded")
			return
		plt.plot(
			*self.interpolateRhoIJ(i1, i2, y, ri=0),
			ls=ls, c=lc, label="%s %d%d re"%(self.label, i1+1, i2+1)
			)
		if im:
			plt.plot(
				*self.interpolateRhoIJ(i1, i2, y, ri=1),
				ls=":", c=lc, label="%s %d%d im"%(self.label, i1+1, i2+1)
				)
		plt.xscale("log")
		plt.xlabel("$x$")
		plt.ylabel(r"$\rho_{ij}$")

	def plotdRhoOffDiagY(self, i1, i2, y, lc="k", ls="-", im=True):
		if not self.full:
			print("no offdiagonal loaded")
			return
		dijrex, dijrey = self.interpolateRhoIJ(i1, i2, y, ri=0)
		plt.plot(
			dijrex, np.gradient(dijrey, dijrex),
			ls=ls, c=lc, label="%s %d%d re"%(self.label, i1+1, i2+1)
			)
		if im:
			dijimx, dijimy = self.interpolateRhoIJ(i1, i2, y, ri=1)
			plt.plot(
				dijimx, np.gradient(dijimy, dijimx),
				ls=":", c=lc, label="%s %d%d im"%(self.label, i1+1, i2+1)
				)
		plt.xscale("log")
		plt.xlabel("$x$")
		plt.ylabel(r"$d\rho_{ij}/dt$")

	def doAllPlots(self, yref=5., color="k"):
		plt.close()
		self.plotZ(lc=color)
		finalizePlot(
			"%s/z.pdf"%self.folder,
			xlab="$x$",
			ylab=r"$z$",
			xscale="log",
			)

		for i in range(self.nnu):
			self.plotRhoDiagY(i, yref, styles[i], lc=colors[i])
		finalizePlot(
			"%s/rho_diag.pdf"%self.folder,
			xlab="$x$",
			ylab=r"$\rho/\rho_{eq}-1$",
			xscale="log",
			)

		for i in range(self.nnu):
			self.plotRhoFin(i, ls=styles[i], lc=colors[i])
		finalizePlot(
			"%s/rhofin_diag.pdf"%self.folder,
			)

		if self.full:
			for i in range(self.nnu):
				for j in range(i+1, self.nnu):
					self.plotRhoOffDiagY(i, j, yref, lc=colors[2*i+j-1])
			finalizePlot(
				"%s/rho_offdiag.pdf"%self.folder,
				)

			for i in range(self.nnu):
				for j in range(i+1, self.nnu):
					self.plotdRhoOffDiagY(i, j, yref, lc=colors[2*i+j-1])
			finalizePlot(
				"%s/drho_offdiag.pdf"%self.folder,
				)
