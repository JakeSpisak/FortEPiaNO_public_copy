"""Functions and classes that read the output and help to do plots"""

import os
import re
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import interp1d
from scipy.integrate import quad
from matplotlib.ticker import AutoMinorLocator

try:
	FileNotFoundError
except NameError:
	FileNotFoundError = IOError

colors = ["r", "g", "b", "k", "c", "m", "y", "#99ff33", "#ff9933"] *4
styles = ["-", "--", ":", "-."] *2
markers = [".", "+", "x", "^", "*", "h", "D"]

PISQD15 = np.pi**2/15.

def finalizePlot(fname,
		lloc="best",
		title="",
		xlab=None,
		ylab=None,
		xscale=None,
		yscale=None,
		xlim=None,
		ylim=None,
		legcol=1,
		legend=True,
		x_T=None,
		Neff_axes=False,
		tightrect=(-0.025, -0.025, 1.02, 1.02),
		):
	plt.title(title)
	ax=plt.gca()
	if not Neff_axes:
		ax.tick_params("both", which="both", direction="in",
			left=True, right=True, top=True, bottom=True)
	if legend:
		plt.legend(loc=lloc, ncol=legcol)
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
	if x_T:
		lims = ax.get_xlim()
		ax.set_xscale("log")
		ax.set_xlabel("$x$")
		ax1=ax.twiny()
		ax1.set_xlim([0.5109989461/lims[0], 0.5109989461/lims[1]])
		ax1.set_xscale("log")
		ax1.set_xlabel("$T$ [MeV]")
	if Neff_axes:
		ax.set_ylabel(r"$N_{\rm eff}^{\rm in}=\frac{8}{7}\frac{\rho_\nu}{\rho_\gamma}$")
		lims = ax.get_ylim()
		ax1=ax.twinx()
		ax.tick_params("both", which="both", direction="out",
			left=True, right=False, labelleft=True, labelright=False)
		ax1.tick_params("both", which="both", direction="out",
			left=False, right=True, labelleft=False, labelright=True)
		ax1.set_ylabel(r"$N_{\rm eff}^{\rm now}=\frac{8}{7}\left(\frac{11}{4}\right)^{4/3}\;\frac{\rho_\nu}{\rho_\gamma}$")
		ax1.set_ylim(np.asarray(lims)*(11./4)**(4./3))
		minorLocatorX = AutoMinorLocator(2)
		ax1.yaxis.set_minor_locator(minorLocatorX)
	plt.tight_layout(rect=tightrect)
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

class FortEPiaNORun():
	"""Class that reads the output and helps to do plots"""

	def __init__(self,
			folder,
			nnu=3,
			full=True,
			label="",
			plots=False,
			rho=True,
			verbose=True,
			):
		self.folder = folder
		self.full = full
		self.label = label
		self.verbose = verbose
		if not os.path.exists(folder):
			if verbose:
				print("non-existing folder: %s"%folder)
			return
		try:
			fdy = np.loadtxt("%s/fd.dat"%folder)
		except (IOError, OSError):
			self.yv = np.nan
			self.fd = np.nan
		else:
			self.yv = fdy[:,0]
			self.fd = fdy[:,1]
		try:
			self.zdat = np.loadtxt("%s/z.dat"%folder)
		except (IOError, OSError):
			self.zdat = np.asarray([[np.nan, np.nan, np.nan]])
		try:
			self.Neffdat = np.loadtxt("%s/Neff.dat"%folder)
		except (IOError, OSError):
			self.Neffdat = np.asarray([[np.nan, np.nan, np.nan]])
		try:
			self.endens = np.loadtxt("%s/energyDensity.dat"%folder)
		except (IOError, OSError):
			self.endens = np.nan
		try:
			self.entropy = np.loadtxt("%s/entropy.dat"%folder)
		except (IOError, OSError):
			self.entropy = np.nan
		self.nnu = nnu
		self.rho = np.asarray([
			[[None, None] for i in range(nnu)] for j in range(nnu)
			])
		self.rhoM = np.asarray([
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
			try:
				self.wfin = float(
					re.match("final w[ =]*([-\d.]*)", self.resume[0]).group(1))
			except AttributeError:
				if verbose:
					print("final w is not in resume.dat")
				zlineindex=0
			else:
				zlineindex = 1
			self.zfin = float(
				re.match("final z[ =]*([-\d.]*)", self.resume[zlineindex]).group(1))
		self.deltarhofin = []
		for i in range(self.nnu):
			if self.hasResume:
				self.deltarhofin.append(
					float(
						re.match("dRho_%s[ =]*([-\d.]*)"%(i+1), self.resume[i+1+zlineindex])
							.group(1)
						)
					)
			try:
				self.rho[i, i, 0] = np.loadtxt(
					"%s/nuDens_diag%d.dat"%(folder, i+1))
			except (IOError, OSError):
				self.rho[i, i, 0] = np.nan
			if full:
				for j in range(i+1, self.nnu):
					try:
						self.rho[i, j, 0] = np.loadtxt(
							"%s/nuDens_nd_%d%d_re.dat"%(folder, i+1, j+1))
					except (IOError, OSError):
						self.rho[i, j, 0] = np.nan
					try:
						self.rho[i, j, 1] = np.loadtxt(
							"%s/nuDens_nd_%d%d_im.dat"%(folder, i+1, j+1))
					except (IOError, OSError):
						self.rho[i, j, 1] = np.nan
			try:
				self.rhoM[i, i, 0] = np.loadtxt(
					"%s/nuDens_mass%d.dat"%(folder, i+1))
			except (IOError, OSError):
				self.rhoM[i, i, 0] = np.nan
		self.printTableLine()
		if rho and plots:
			self.doAllPlots()

	def interpolateRhoIJ(self, i1, i2, y, ri=0, y2=False, mass=False):
		if mass:
			rho = self.rhoM
		else:
			rho = self.rho
		xv = []
		yv = []
		prevy = 0
		for i, x in enumerate(rho[i1, i2, ri][:, 0]):
			fy = interp1d(self.yv, rho[i1, i2, ri][i, 1:])
			cy = fy(y)
			if cy != prevy:
				prevy = cy
				yv.append(prevy * (y**2 if y2 else 1.))
				xv.append(x)
		xv.append(x)
		yv.append(cy)
		return xv, yv

	def interpolateRhoIJ_x(self, i1, i2, x, ri=0, y2=False, mass=False):
		if mass:
			rho = self.rhoM
		else:
			rho = self.rho
		ov = []
		for i, y in enumerate(self.yv):
			fx = interp1d(
				rho[i1, i2, ri][:, 0],
				rho[i1, i2, ri][:, i] * (y**2 if y2 else 1.))
			ov.append(fx(x))
		return self.yv, ov

	def printTableLine(self):
		if not self.verbose:
			return
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
		else:
			print("{lab:<35s} \t\t\t\tnot finished, currently on x={x:}".format(
				lab=self.label,
				x=self.zdat[-1],
				))

	def plotFD(self, ls="-", lc="k", lab=None, rescale=1., fac=1.):
		if rescale != 1.:
			fd = self.fd*(np.exp(self.yv)+1.)/(np.exp(self.yv/rescale)+1.)
		else:
			fd = self.fd
		plt.plot(
			self.yv, fac*fd,
			label=self.label if lab is None else lab, ls=ls, marker=".", c=lc
			)
		plt.xscale("log")
		plt.yscale("log")
		plt.xlabel("$y$")
		plt.ylabel(r"$y^2 f(y)$")

	def plotZ(self, ls="-", lc="k", lab=None):
		plt.plot(
			*stripRepeated(self.zdat, 0, 1),
			label=self.label if lab is None else lab, ls=ls, c=lc
			)
		plt.xscale("log")
		plt.xlabel("$x$")
		plt.ylabel(r"$z$")

	def plotW(self, ls="-", lc="k", lab=None):
		try:
			self.zdat[0, 2]
		except IndexError:
			print("w is not in z.dat")
			return
		plt.plot(
			*stripRepeated(self.zdat, 0, 2),
			label=self.label if lab is None else lab, ls=ls, c=lc
			)
		plt.xscale("log")
		plt.xlabel("$x$")
		plt.ylabel(r"$w$")

	def plotZoverW(self, ls="-", lc="k", lab=None):
		try:
			self.zdat[0, 2]
		except IndexError:
			print("w is not in z.dat")
			return
		plt.plot(
			*stripRepeated(np.asarray([[x[0], x[1]/x[2]] for x in self.zdat]), 0, 1),
			label=self.label if lab is None else lab, ls=ls, c=lc
			)
		plt.xscale("log")
		plt.xlabel("$x$")
		plt.ylabel(r"$w$")

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

	def plotRhoDiag(self, inu, iy, ls, lc="k", mass=False):
		if mass:
			rho = self.rhoM
		else:
			rho = self.rho
		plt.plot(
			*stripRepeated(rho[inu, inu, 0], 0, iy),
			label="%s \alpha=%d"%(self.label, inu+1), ls=ls, c=lc
			)
		plt.xscale("log")
		plt.xlabel("$x$")
		plt.ylabel(r"$\rho_{\alpha\alpha}$")

	def plotdRhoDiag(self, inu, iy, ls, lc="k", mass=False):
		if mass:
			rho = self.rhoM
		else:
			rho = self.rho
		dijrex, dijrey = stripRepeated(rho[inu, inu, 0], 0, iy)
		plt.plot(
			dijrex, np.gradient(dijrey, dijrex),
			ls=ls, c=lc, label="%s \alpha=%d"%(self.label, inu+1)
			)
		plt.xscale("log")
		plt.xlabel("$x$")
		plt.ylabel(r"$d\rho_{\alpha\alpha}/dt$")

	def plotRhoOffDiag(self, i1, i2, iy, lc="k", im=True, mass=False):
		if not self.full:
			print("no offdiagonal loaded")
			return
		if mass:
			rho = self.rhoM
		else:
			rho = self.rho
		plt.plot(
			*stripRepeated(rho[i1, i2, 0], 0, iy),
			ls="-", c=lc, label="%s \alpha\beta=%d%d re"%(self.label, i1+1, i2+1)
			)
		if im:
			plt.plot(
				*stripRepeated(rho[i1, i2, 1], 0, iy),
				ls=":", c=lc, label="%s \alpha\beta=%d%d im"%(self.label, i1+1, i2+1)
				)
		plt.xscale("log")
		plt.xlabel("$x$")
		plt.ylabel(r"$\rho_{\alpha\beta}$")

	def plotdRhoOffDiag(self, i1, i2, iy, lc="k", im=True, mass=False):
		if not self.full:
			print("no offdiagonal loaded")
			return
		if mass:
			rho = self.rhoM
		else:
			rho = self.rho
		dijrex, dijrey = stripRepeated(rho[i1, i2, 0], 0, iy)
		plt.plot(
			dijrex, np.gradient(dijrey, dijrex),
			ls="-", c=lc, label="%s \alpha\beta=%d%d re"%(self.label, i1+1, i2+1)
			)
		if im:
			dijimx, dijimy = stripRepeated(rho[i1, i2, 1], 0, iy)
			plt.plot(
				dijimx, np.gradient(dijimy, dijimx),
				ls=":", c=lc, label="%s \alpha\beta=%d%d im"%(self.label, i1+1, i2+1)
				)
		plt.xscale("log")
		plt.xlabel("$x$")
		plt.ylabel(r"$d\rho_{\alpha\beta}/dt$")

	def plotRhoFin(self, ix, iy=None, ri=0, ls="-", lc="k", y2=False, lab=None, mass=False):
		if mass:
			rho = self.rhoM
		else:
			rho = self.rho
		if iy is None:
			iy = ix
		if ri not in [0, 1]:
			ri = 0
		label = "%s \alpha\beta=%d%d %s"%(self.label, ix+1, iy+1, "re" if ri == 0 else "im") \
			if lab is None else lab
		fyv = self.yv**2*rho[ix, iy, ri][-1, 1:] if y2 else rho[ix, iy, ri][-1, 1:]
		plt.plot(
			self.yv, fyv,
			ls=ls, c=lc,
			label=label,
			)
		plt.xlabel("$y$")
		plt.ylabel(r"$%s\rho_{\alpha\beta}^{\rm fin}(y)$"%("y^2" if y2 else ""))

	def plotRhoX(self, i1, x, i2=None, ri=0, ls="-", lc="k", y2=False, mass=False):
		if i2 is None:
			i2 = i1
		if ri not in [0, 1]:
			ri = 0
		plt.plot(
			*self.interpolateRhoIJ_x(i1, i2, x, ri, y2=y2, mass=mass),
			ls=ls, c=lc,
			label="%s \alpha\beta=%d%d %s x=%f"%(self.label, i1+1, i2+1, "re" if ri == 0 else "im", x)
			)
		plt.xlabel("$y$")
		plt.ylabel(r"$%s\rho_{\alpha\beta}(y)$"%("y^2" if y2 else ""))

	def plotRhoDiagY(self, inu, y, ls, lc="k", lab=None, y2=False, mass=False):
		x, yv = self.interpolateRhoIJ(inu, inu, y, ri=0, mass=mass)
		label = lab if lab is not None else "%s \alpha=%d"%(self.label, inu+1)
		plt.plot(
			x, np.asarray(yv) * (y**2 if y2 else 1.),
			label=label, ls=ls, c=lc
			)
		plt.xscale("log")
		plt.xlabel("$x$")
		plt.ylabel(r"$%s\rho_{\alpha\alpha}$"%("y^2" if y2 else ""))

	def plotdRhoDiagY(self, inu, y, ls, lc="k", lab=None, y2=False, mass=False):
		x, yv = self.interpolateRhoIJ(inu, inu, y, ri=0, mass=mass)
		label = lab if lab is not None else "%s \alpha=%d"%(self.label, inu+1)
		plt.plot(
			x, np.gradient(np.asarray(yv) * (y**2 if y2 else 1.), x),
			ls=ls, c=lc, label="%s \alpha=%d"%(self.label, inu+1)
			)
		plt.xscale("log")
		plt.xlabel("$x$")
		plt.ylabel(r"$d\rho_{\alpha\alpha}/dt$")

	def plotRhoOffDiagY(self, i1, i2, y, lc="k", ls="-", im=True, lab=None, mass=False):
		if not self.full:
			print("no offdiagonal loaded")
			return
		plt.plot(
			*self.interpolateRhoIJ(i1, i2, y, ri=0, mass=mass),
			ls=ls, c=lc, label="%s \alpha\beta=%d%d re"%(self.label, i1+1, i2+1) if lab is None else lab
			)
		if im:
			plt.plot(
				*self.interpolateRhoIJ(i1, i2, y, ri=1),
				ls=":", c=lc, label="%s \alpha\beta=%d%d im"%(self.label, i1+1, i2+1) if lab is None else lab
				)
		plt.xscale("log")
		plt.xlabel("$x$")
		plt.ylabel(r"$\rho_{\alpha\beta}$")

	def plotdRhoOffDiagY(self, i1, i2, y, lc="k", ls="-", im=True, lab=None, mass=False):
		if not self.full:
			print("no offdiagonal loaded")
			return
		dijrex, dijrey = self.interpolateRhoIJ(i1, i2, y, ri=0, mass=mass)
		try:
			plt.plot(
				dijrex, np.gradient(dijrey, dijrex),
				ls=ls, c=lc, label="%s %d%d re"%(self.label, i1+1, i2+1) if lab is None else lab
				)
		except IndexError:
			pass
		if im:
			dijimx, dijimy = self.interpolateRhoIJ(i1, i2, y, ri=1, mass=mass)
			try:
				plt.plot(
					dijimx, np.gradient(dijimy, dijimx),
					ls=":", c=lc, label="%s %d%d im"%(self.label, i1+1, i2+1) if lab is None else lab
					)
			except IndexError:
				pass
		plt.xscale("log")
		plt.xlabel("$x$")
		plt.ylabel(r"$d\rho_{\alpha\beta}/dt$")

	def plotNeff(self, lc="k", ls="-", im=True, lab=None, xlims=[0.5, 4.5], axes=True):
		if not np.isnan(self.Neffdat[0,0]):
			data = self.Neffdat
		else:
			data = []
			for ix, [x, z] in enumerate(self.zdat[:, 0:2]):
				rhogamma = PISQD15 * z**4
				rhonu = np.sum([self.integrateRho_yn(inu, 3, ix=ix) for inu in range(self.nnu)])
				data.append([x, 8./7.*rhonu/rhogamma, 8./7.*rhonu/rhogamma*(11./4.)**(4./3.)])
			data = np.asarray(data)
			print(os.path.join(self.folder, "Neff.dat"))
			np.savetxt(os.path.join(self.folder, "Neff.dat"), data, fmt="%.7e")
		plt.plot(
			*stripRepeated(data, 0, 1),
			ls=ls, c=lc, label=self.label if lab is None else lab
			)
		plt.xscale("log")
		plt.xlabel("$x$")
		if axes:
			ax=plt.gca()
			ax.set_ylim(xlims)
			lims = ax.get_ylim()
			ax1=ax.twinx()
			ax.set_ylabel(r"$N_{\rm eff}^{\rm in}$")
			ax1.set_ylabel(r"$N_{\rm eff}^{\rm now}$")
			ax1.set_ylim(np.asarray(lims)*(11./4)**(4./3))

	def plotEnergyDensity(self,
			gamma_e=True,
			gec="#00ccff",
			ges=":",
			gamma_e_mu=True,
			gemc="#6666ff",
			gems="--",
			labels=[r"$\gamma$", "$e$", r"$\mu$", r"$\nu_e$", r"$\nu_\mu$", r"$\nu_\tau$", r"$\nu_s$"],
			colors=["r", "b", "g", "#ff9933", "#ff9933", "#ff9933", "#ff00ff"],
			styles=["-", "-", "-", ":", "-.", "--", "-"],
			skip=[False, False, False, False, False, False, False],
			):
		plt.plot(self.endens[:,0], np.asarray([np.sum(cl[2:]) for cl in self.endens]), label="total", c="k")
		for ix, lab in enumerate(labels):
			if skip[ix]:
				continue
			try:
				plt.plot(self.endens[:, 0], self.endens[:, 2+ix], label=lab, c=colors[ix], ls=styles[ix])
			except IndexError:
				pass
		if gamma_e:
			plt.plot(
				self.endens[:,0],
				self.endens[:,2]+self.endens[:,3],
				label=r"$\gamma+e$",
				c=gec,
				ls=ges,
				)
		if gamma_e_mu:
			plt.plot(
				self.endens[:,0],
				self.endens[:,2]+self.endens[:,3]+self.endens[:,4],
				label=r"$\gamma+e+\mu$",
				c=gemc,
				ls=gems,
				)

	def plotEntropy(self,
			gamma_e=True,
			gec="#00ccff",
			ges=":",
			gamma_e_mu=True,
			gemc="#6666ff",
			gems="--",
			labels=[r"$\gamma$", "$e$", r"$\mu$", r"$\nu_e$", r"$\nu_\mu$", r"$\nu_\tau$", r"$\nu_s$"],
			colors=["r", "b", "g", "#ff9933", "#ff9933", "#ff9933", "#ff00ff"],
			styles=["-", "-", "-", ":", "-.", "--", "-"],
			skip=[False, False, False, False, False, False, False],
			lw=1,
			allstyles=False,
			alllabels=None,
			):
		plt.plot(
			self.entropy[:,0],
			np.asarray([np.sum(cl[2:]) for cl in self.entropy]),
			label="total" if alllabels is None else alllabels,
			c="k",
			ls="-" if not allstyles else allstyles,
			lw=lw,
			)
		for ix, lab in enumerate(labels):
			if skip[ix]:
				continue
			try:
				plt.plot(
					self.entropy[:, 0],
					self.entropy[:, 2+ix],
					label=lab if alllabels is None else alllabels,
					c=colors[ix],
					ls=styles[ix] if not allstyles else allstyles,
					lw=lw,
					)
			except IndexError:
				pass
		if gamma_e:
			plt.plot(
				self.entropy[:,0],
				self.entropy[:,2]+self.entropy[:,3],
				label=r"$\gamma+e$" if alllabels is None else alllabels,
				c=gec,
				ls=ges if not allstyles else allstyles,
				lw=lw,
				)
		if gamma_e_mu:
			plt.plot(
				self.entropy[:,0],
				self.entropy[:,2]+self.entropy[:,3]+self.entropy[:,4],
				label=r"$\gamma+e+\mu$" if alllabels is None else alllabels,
				c=gemc,
				ls=gems if not allstyles else allstyles,
				lw=lw,
				)

	def doAllPlots(self, yref=5., color="k"):
		plt.close()
		self.plotZ(lc=color, lab="z")
		self.plotW(lc=color, ls=":", lab="w")
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
			ylab=r"$\rho$",
			xscale="log",
			yscale="log",
			)

		for i in range(self.nnu):
			self.plotRhoDiagY(i, yref, styles[i], lc=colors[i], mass=True)
		finalizePlot(
			"%s/rho_mass_diag.pdf"%self.folder,
			xlab="$x$",
			ylab=r"$\rho$",
			xscale="log",
			yscale="log",
			)

		for i in range(self.nnu):
			self.plotRhoFin(i, ls=styles[i], lc=colors[i], y2=True)
		finalizePlot(
			"%s/rhofin_diag.pdf"%self.folder,
			xscale="linear",
			yscale="log",
			)

		for i in range(self.nnu):
			self.plotRhoFin(i, ls=styles[i], lc=colors[i], y2=True, mass=True)
		finalizePlot(
			"%s/rhofin_mass_diag.pdf"%self.folder,
			xscale="linear",
			yscale="log",
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

	def integrateRho_yn(self, inu, n, ix=-1, show=False, mass=False):
		"""Compute the integral
		Int_0^Inf dy y^n f(y)/Pi^2
		for the requested eigenstate at the given x
		"""
		if mass:
			rho = self.rhoM
		else:
			rho = self.rho
		fy = interp1d(self.yv, rho[inu, inu, 0][ix, 1:]*(np.exp(self.yv)+1))
		res = quad(lambda y: y**n * fy(y)/(np.exp(y)+1), 0.01, 20)
		if show:
			print(res)
		return res[0]/np.pi**2
