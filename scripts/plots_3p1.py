import traceback
import matplotlib
matplotlib.use("agg")
import matplotlib.pyplot as plt
import numpy as np
from fortepianoOutput import colors, styles, finalizePlot, stripRepeated, FortEPiaNORun

flavors=["e", r"\mu", r"\tau", "s"]
a3nu = FortEPiaNORun("OUT/3nu/n30_hp/", label=r"$3\nu$")
# new3p1 = FortEPiaNORun("OUT/prec_3p1/n40_hp/", nnu=4, label=r"3+1")
Ue4h = FortEPiaNORun("OUT/3p1/30_1_0.01_0_0", nnu=4, label=r"$|U_{e4}|^2=10^{-2}$")
Ue4l = FortEPiaNORun("OUT/3p1/30_1_0.001_0_0", nnu=4, label=r"$|U_{e4}|^2=10^{-3}$")
Um4h = FortEPiaNORun("OUT/3p1/30_1_0_0.01_0", nnu=4, label=r"$|U_{\mu4}|^2=10^{-2}$")
Um4l = FortEPiaNORun("OUT/3p1/30_1_0_0.0001_0", nnu=4, label=r"$|U_{\mu4}|^2=10^{-4}$")
Ut4h = FortEPiaNORun("OUT/3p1/30_1_0_0_0.01", nnu=4, label=r"$|U_{\tau4}|^2=10^{-2}$")
Ut4l = FortEPiaNORun("OUT/3p1/30_1_0_0_0.0001", nnu=4, label=r"$|U_{\tau4}|^2=10^{-4}$")

# plot Neff with three different choices for the angles
fig = plt.figure(figsize=(6,4))
for ir, r in enumerate([
		Ue4h, Um4h, Ut4h, None, Ue4l, a3nu]):
	if r is not None:
		print("Neff for %s"%r.label)
		r.plotNeff(lc=colors[ir], ls='-', axes=False)
finalizePlot(
	"plots/3p1/Neff.pdf",
	xlim=[1e-3, 35],
	ylim=[0.5, 4.5],
	xlab="$x$",
	legend=True,
	lloc="lower left",
	x_T=True,
	Neff_axes=True,
	)

#energy densities
fig = plt.figure(figsize=(5.1,3.4))
a3nu.plotEnergyDensity(
	labels=[r"$\gamma$", "$e$", r"$\mu$", r"$\nu_e,\nu_\mu,\nu_\tau$", "", "", r"$\nu_s$"],
	skip=[False, False, False, False, True, True, False],
	allstyles=":",
	alllabels="",
	lw=1,
	)
Um4h.plotEnergyDensity(
	labels=[r"$\gamma$", "$e$", r"$\mu$", r"$\nu_e,\nu_\mu,\nu_\tau$", "", "", r"$\nu_s$"],
	skip=[False, False, False, False, True, True, False],
	allstyles="-",
	alllabels=None,
	lw=2,
	)
finalizePlot(
	"plots/3p1/energyDensity.pdf",
	xlab="$x$",
	xlim=[1e-3, 35],
	ylim=[0, 8.5],
	ylab=r"$\rho$",
	xscale="log",
	x_T=True,
	legcol=4,
	)

#entropy combined
fig = plt.figure(figsize=(5.1,3.4))
a3nu.plotEntropy(
	labels=[r"$\gamma$", "$e$", r"$\mu$", r"$\nu_e,\nu_\mu,\nu_\tau$", "", "", r"$\nu_s$"],
	skip=[False, False, False, False, True, True, False],
	allstyles=":",
	alllabels="",
	lw=1,
	)
Um4h.plotEntropy(
	labels=[r"$\gamma$", "$e$", r"$\mu$", r"$\nu_e$, $\nu_\mu$, $\nu_\tau$", "", "", r"$\nu_s$"],
	skip=[False, False, False, False, True, True, False],
	allstyles="-",
	alllabels=None,
	lw=2,
	)
finalizePlot(
	"plots/3p1/entropy_combined.pdf",
	xlab="$x$",
	xlim=[1e-3, 35],
	ylim=[0, 9],
	ylab=r"$s$",
	xscale="log",
	x_T=True,
	legcol=4,
	)

# rho diagonal for different y, only one run
for iy, y in enumerate([0.1, 2, 5]):
	plt.plot(np.nan, ls='-', color=colors[iy], label="y=%s"%y)
for i in [0, 3]:
	plt.plot(np.nan, ls=styles[i], color="k", label=r"$\alpha=%s$"%flavors[i])
for iy, y in enumerate([0.3, 2, 5]):
	for i in [0, 3]:
		Ue4h.plotRhoDiagY(i, y, styles[i], lc=colors[iy], lab="", y2=True)
finalizePlot(
	"plots/3p1/rho_diag_y.pdf",
	xlab="$x$",
	ylab=r"$y^2\rho_{\alpha\alpha}$",
	xlim=(1e-3, 30),
	ylim=(1e-5, 1.),
	xscale="log",
	yscale="log",
	lloc="lower right",
	legcol=2,
	x_T=True,
	)

# evolution of rho(y=5) for different cases
for i in range(4):
	plt.plot(np.nan, ls=styles[i], color="k", label=r"$\alpha=%s$"%flavors[i])
plt.plot(np.nan, ls='-', color=colors[0], label=Ue4h.label)
plt.plot(np.nan, ls='-', color=colors[1], label=Ue4l.label)
plt.plot(np.nan, ls='-', color=colors[2], label=a3nu.label)
for i in range(4):
	Ue4h.plotRhoDiagY(i, 5, styles[i], lc=colors[0], y2=True, lab="")
	Ue4l.plotRhoDiagY(i, 5, styles[i], lc=colors[1], y2=True, lab="")
for i in range(3):
	a3nu.plotRhoDiagY(i, 5, styles[i], lc=colors[2], y2=True, lab="")
finalizePlot(
	"plots/3p1/rho_diag_prec.pdf",
	xlab="$x$",
	ylab=r"$y^2\rho_{\alpha\alpha}$",
	xlim=(1e-3, 30),
	ylim=(1e-3, 0.3),
	xscale="log",
	yscale="log",
	lloc="lower right",
	x_T=True,
	legcol=2,
	)

# plots of z and w for some cases
for ir, r in enumerate([Ue4h, Ue4l, a3nu]):
	plt.plot(np.nan, ls='-', color=colors[ir], label=r.label)
	r.plotZ(lc=colors[ir], ls='-', lab="")
	r.plotW(lc=colors[ir], ls=':', lab="")
plt.plot(np.nan, ls='-', color=colors[3], label=r"$z$")
plt.plot(np.nan, ls=':', color=colors[3], label=r"$w$")
finalizePlot(
	"plots/3p1/z_prec.pdf",
	xlim=[1e-3, 35],
	xlab="$x$",
	ylab=r"$w,z$",
	xscale="log",
	x_T=True,
	legcol=2,
	)

# offdiagonal rho and drho/dx for two cases
for ol, fn in [
		[[a3nu, Ue4h], "prec"],
		]:
	for i in range(4):
		for j in range(i+1, 4):
			plt.plot(np.nan, ls='-', color=colors[2*i+j-1], label=r"$\alpha\beta=%s %s$"%(flavors[i], flavors[j]))
	first_legend = plt.legend(loc='lower left')
	ax = plt.gca()
	ax.add_artist(first_legend)
	for io, o in enumerate(ol):
		plt.plot(np.nan, ls=styles[io], color='k', label=o.label)
	hands, labs = ax.get_legend_handles_labels()
	second_legend = plt.legend(hands[-2:], labs[-2:], loc='lower right')
	for io, o in enumerate(ol):
		for i in range(4):
			for j in range(i+1, 4):
				try:
					o.plotRhoOffDiagY(i, j, 5, ls=styles[io], lc=colors[2*i+j-1], im=False, lab="")
				except IndexError:
					pass
	finalizePlot(
		"plots/3p1/rho_offdiag_%s.pdf"%fn,
		legend=False,
		x_T=True,
		)

	for i in range(4):
		for j in range(i+1, 4):
			plt.plot(np.nan, ls='-', color=colors[2*i+j-1], label=r"$\alpha\beta=%s %s$"%(flavors[i], flavors[j]))
	first_legend = plt.legend(loc='lower left')
	ax = plt.gca()
	ax.add_artist(first_legend)
	for io, o in enumerate(ol):
		plt.plot(np.nan, ls=styles[io], color='k', label=o.label)
	hands, labs = ax.get_legend_handles_labels()
	second_legend = plt.legend(hands[-2:], labs[-2:], loc='lower right')
	for io, o in enumerate(ol):
		for i in range(4):
			for j in range(i+1, 4):
				try:
					o.plotdRhoOffDiagY(i, j, 5, ls=styles[io], lc=colors[2*i+j-1], im=False, lab="")
				except IndexError:
					pass
	finalizePlot(
		"plots/3p1/rho_doffdiag_%s.pdf"%fn,
		legend=False,
		x_T=True,
		)

# final rho for some cases
for i in range(4):
	plt.plot(np.nan, ls=styles[i], color="k", label=r"$\alpha=%s$"%flavors[i])
for fn, cl, ol in [
		["ord", (0, 1, 2), (Ue4h, Ue4l, a3nu)],
		]:
	for o, c in zip(ol, cl):
		plt.plot(np.nan, ls='-', color=colors[c], label=o.label)
		for i in range(4):
			try:
				o.plotRhoFin(i, ls=styles[i], lc=colors[c], y2=True, lab="")
			except IndexError:
				pass#print(traceback.format_exc())
	a3nu.plotFD(lab=r"FD $3+1$ DW", ls=":", lc=colors[3], rescale=Ue4l.wfin, fac=(Ue4l.Neff-a3nu.Neff))
	finalizePlot(
		"plots/3p1/rhofin_diag_%s.pdf"%fn,
		lloc="upper right",
		xlim=[-0.2, 15],
		ylim=[1e-4, 1],
		ylab=r"$y^2\rho_{\alpha\alpha}^{\rm fin}(y)$",
		xscale="linear",
		yscale="log",
		legcol=2,
		)

# plot z and w with three different choices for the angles
for ir, r in enumerate([Ue4h, Um4h, Ut4h, None, Ue4l, a3nu]):
	if r is not None:
		plt.plot(np.nan, ls='-', color=colors[ir], label=r.label)
		r.plotZ(lc=colors[ir], ls='-', lab="")
		r.plotW(lc=colors[ir], ls=':', lab="")
first_legend = plt.legend(loc='upper left')
ax = plt.gca()
ax.add_artist(first_legend)
plt.plot(np.nan, ls='-', color=colors[3], label=r"$z$")
plt.plot(np.nan, ls=':', color=colors[3], label=r"$w$")
a3nu.plotZoverW(ls="--", lc="#999999", lab="$z/w$")
hands, labs = ax.get_legend_handles_labels()
second_legend = plt.legend(hands[-3:], labs[-3:], loc='upper center')
finalizePlot(
	"plots/3p1/z_angles.pdf",
	xlim=[1e-3, 35],
	xlab="$x$",
	ylab=r"$w,z$",
	xscale="log",
	legend=False,
	x_T=True,
	)

# evolution of rho(y=5) for different cases
for i in range(4):
	plt.plot(np.nan, ls=styles[i], color="k", label=r"$\alpha=%s$"%flavors[i])
first_legend = plt.legend(loc='lower left')
ax = plt.gca()
ax.add_artist(first_legend)
for ir, r in enumerate([Ue4l, Um4l, Ue4l, None, Ue4h, Um4h, Ut4h]):
	if r is not None:
		plt.plot(np.nan, ls='-', color=colors[ir], label=r.label)
		for i in range(4):
			r.plotRhoDiagY(i, 5, styles[i], lc=colors[ir], y2=True, lab="")
hands, labs = ax.get_legend_handles_labels()
second_legend = plt.legend(hands[4:], labs[4:], loc='lower right', ncol=2)
finalizePlot(
	"plots/3p1/rho_diag_angles.pdf",
	xlab="$x$",
	ylab=r"$y^2\rho_{\alpha\alpha}$",
	xlim=(1e-3, 30),
	ylim=(1e-3, 0.3),
	xscale="log",
	yscale="log",
	legend=False,
	x_T=True,
	)

# # offdiagonal rho and drho/dx for two cases
for ol, fn in [
		[[Ue4h, Um4h, Ut4h], "angles"],
		]:
	for i in range(4):
		for j in range(i+1, 4):
			plt.plot(np.nan, ls='-', color=colors[2*i+j-1], label=r"$\alpha\beta=%s %s$"%(flavors[i], flavors[j]))
	first_legend = plt.legend(loc='lower left')
	ax = plt.gca()
	ax.add_artist(first_legend)
	for io, o in enumerate(ol):
		plt.plot(np.nan, ls=styles[io], color='k', label=o.label)
	hands, labs = ax.get_legend_handles_labels()
	second_legend = plt.legend(hands[-3:], labs[-3:], loc='lower right')
	for io, o in enumerate(ol):
		for i in range(4):
			for j in range(i+1, 4):
				try:
					o.plotRhoOffDiagY(i, j, 5, ls=styles[io], lc=colors[2*i+j-1], im=False, lab="")
				except IndexError:
					pass
	finalizePlot(
		"plots/3p1/rho_offdiag_%s.pdf"%fn,
		xlim=(1e-3, 30),
		legend=False,
		lloc="lower right", legcol=3,
		x_T=True,
		)

print("IO:")
Ue4ih = FortEPiaNORun("OUT/3p1/30i_1_0.01_0_0", nnu=4, label=r"$|U_{e4}|^2=10^{-2}$")
Ue4il = FortEPiaNORun("OUT/3p1/30i_1_0.001_0_0", nnu=4, label=r"$|U_{e4}|^2=10^{-3}$")
Um4ih = FortEPiaNORun("OUT/3p1/30i_1_0_0.01_0", nnu=4, label=r"$|U_{\mu4}|^2=10^{-2}$")
Um4il = FortEPiaNORun("OUT/3p1/30i_1_0_0.0001_0", nnu=4, label=r"$|U_{\mu4}|^2=10^{-4}$")
Ut4ih = FortEPiaNORun("OUT/3p1/30i_1_0_0_0.01", nnu=4, label=r"$|U_{\tau4}|^2=10^{-2}$")
Ut4il = FortEPiaNORun("OUT/3p1/30i_1_0_0_0.0001", nnu=4, label=r"$|U_{\tau4}|^2=10^{-4}$")

# plot Neff with three different choices for the angles
fig = plt.figure(figsize=(6,4))
for ir, r in enumerate([
		Ue4h, Um4h, Ut4h, Ue4ih, Um4ih, Ut4ih]):
	if r is not None:
		print("Neff for %s"%r.label)
		r.plotNeff(lc=colors[ir], ls='-', axes=False)
finalizePlot(
	"plots/3p1/Neff_ord.pdf",
	xlim=[1e-3, 35],
	ylim=[0.5, 4.5],
	xlab="$x$",
	legend=True,
	lloc="lower left",
	x_T=True,
	Neff_axes=True,
	)
