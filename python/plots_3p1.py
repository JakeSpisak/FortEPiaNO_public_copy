import traceback
import matplotlib
matplotlib.use("agg")
import matplotlib.pyplot as plt
import numpy as np
from nuDensOutput import colors, styles, finalizePlot, stripRepeated, NuDensRun

a3nu = NuDensRun("grids/full_no/OUT/3nu/", label="3nu")
as3p1 = NuDensRun("OUT/prec_3p1/30_hp_4_0.0005/", nnu=4, label=r"3+1")
# grids/full_no/OUT/1.0_0.001_0.001_0.001/
# grids/full_no/OUT/1.0_0.0001_0.0001_0.0001/

for iy, y in enumerate([1, 5, 10]):
	plt.plot([1e3, 2e3], [-10, -11], ls='-', color=colors[iy], label="y=%s"%y)
for i in [0, 3]:
	plt.plot([1e3, 2e3], [-10, -11], ls=styles[i], color="k", label="i=%d"%i)
for iy, y in enumerate([1, 5, 10]):
	for i in [0, 3]:
		as3p1.plotRhoDiagY(i, y, styles[i], lc=colors[iy], p1=1, label="")
finalizePlot(
	"plots/3p1/rho_diag_y.pdf",
	xlab="$x$",
	ylab=r"$\rho_{ii}/\rho_{eq}$",
	xlim=(5e-4, 30),
	ylim=(-0.05, 1.05),
	xscale="log",
	lloc="lower right", legcol=2,
	)

for i in range(3):
	a3nu.plotRhoDiagY(i, 5, styles[i], lc=colors[3],p1=1)
for i in range(4):
	as3p1.plotRhoDiagY(i, 5, styles[i], lc=colors[0], p1=1)
finalizePlot(
	"plots/3p1/rho_diag_prec.pdf",
	xlab="$x$",
	ylab=r"$\rho_{ii}/\rho_{eq}$",
	xscale="log",
	lloc="lower right",
	)

a3nu.plotZ(lc=colors[3])
as3p1.plotZ(lc=colors[0])
finalizePlot(
	"plots/3p1/z_prec.pdf",
	xlim=[0.0005, 35],
	xlab="$x$",
	ylab=r"$z$",
	xscale="log",
	)

for ol, fn in [
		[[a3nu, as3p1], "prec"],
		]:
	for io, o in enumerate(ol):
		for i in range(4):
			for j in range(i+1, 4):
				try:
					o.plotRhoOffDiagY(i, j, 5, ls=styles[io], lc=colors[2*i+j-1], im=False)
				except IndexError:
					pass
	finalizePlot(
		"plots/3p1/rho_offdiag_%s.pdf"%fn,
		lloc="lower right",
		)

	for io, o in enumerate(ol):
		for i in range(4):
			for j in range(i+1, 4):
				try:
					o.plotdRhoOffDiagY(i, j, 5, ls=styles[io], lc=colors[2*i+j-1], im=False)
				except IndexError:
					pass
	finalizePlot(
		"plots/3p1/rho_doffdiag_%s.pdf"%fn,
		lloc="lower right",
		)

for fn, cl, ol in [
		["ord", (3, 0), (a3nu, as3p1)],
		]:
	for o, c in zip(ol, cl):
		for i in range(4):
			try:
				o.plotRhoFin(i, ls=styles[i], lc=colors[c], p1=1)
			except IndexError:
				pass#print(traceback.format_exc())
	finalizePlot(
		"plots/3p1/rhofin_diag_%s.pdf"%fn,
		lloc="lower left",
		# xlim=[0, 12],
		# ylim=[-0.01, 0.07],
		)
