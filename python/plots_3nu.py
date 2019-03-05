import matplotlib
matplotlib.use("agg")
import matplotlib.pyplot as plt
import numpy as np
from nuDensOutput import colors, styles, finalizePlot, stripRepeated, NuDensRun

fNO = NuDensRun("OUT/3nu/20_ci_no/", label="NO")
fIO = NuDensRun("OUT/3nu/20_ci_io/", label="IO")
fNOd = NuDensRun("OUT/3nu/20_df_no/", label="NO damp")
fIOd = NuDensRun("OUT/3nu/20_df_io/", label="IO damp")
fnoosc20 = NuDensRun("OUT/3nu/20_noosc/", label="no osc Ny=20", full=False)
fnoosc40 = NuDensRun("OUT/3nu/40_noosc/", label="no osc Ny=40", full=False)
fnoosc100 = NuDensRun("OUT/3nu/100_noosc/", label="no osc Ny=100", full=False)

for i in range(3):
	fnoosc20.plotRhoDiag(i, 3, styles[i], lc=colors[3])
	fnoosc40.plotRhoDiag(i, 8, styles[i], lc=colors[5])
	fnoosc100.plotRhoDiag(i, 19, styles[i], lc=colors[6])
finalizePlot(
	"plots/3nu/rho_diag_noosc.pdf",
	xlab="$x$",
	ylab=r"$\rho/\rho_{eq}-1$",
	xscale="log",
	lloc="upper left",
	)

fnoosc20.label = "no osc"
for i in range(3):
    fnoosc20.plotRhoDiag(i, 5, styles[i], lc=colors[3])
    fNO.plotRhoDiag(i, 5, styles[i], lc=colors[0])
    fIO.plotRhoDiag(i, 5, styles[i], lc=colors[1])
finalizePlot(
	"plots/3nu/rho_diag_20_ord.pdf",
	xlab="$x$",
	ylab=r"$\rho/\rho_{eq}-1$",
	xscale="log",
	lloc="upper left",
	)

for i in range(3):
	fNO.plotRhoDiag(i, 5, styles[i], lc=colors[0])
	fIO.plotRhoDiag(i, 5, styles[i], lc=colors[1])
	fNOd.plotRhoDiag(i, 5, styles[i], lc=colors[2])
	fIOd.plotRhoDiag(i, 5, styles[i], lc=colors[4])
finalizePlot(
	"plots/3nu/rho_diag_20_damp.pdf",
	xlab="$x$",
	ylab=r"$\rho/\rho_{eq}-1$",
	xscale="log",
	lloc="upper left",
	)

fnoosc20.plotZ(lc=colors[3])
fNO.plotZ(lc=colors[0])
fIO.plotZ(lc=colors[1])
fNOd.plotZ(lc=colors[2])
fIOd.plotZ(lc=colors[4])
finalizePlot(
	"plots/3nu/z_20.pdf",
	xlab="$x$",
	ylab=r"$z$",
	xscale="log",
	)

fNO.plotDeltaZ(fnoosc20, lc=colors[0])
fIO.plotDeltaZ(fnoosc20, lc=colors[1])
fNOd.plotDeltaZ(fnoosc20, lc=colors[2])
fIOd.plotDeltaZ(fnoosc20, lc=colors[4])
finalizePlot(
	"plots/3nu/deltaz_20.pdf",
	xlab="$x$",
	ylab=r"$z-z_{\rm no osc}$",
	xscale="log",
	)

for o1, o2, fn in [
		[fNO, fIO, "ord"],
		[fNO, fNOd, "no"],
		[fIO, fIOd, "io"],
		]:
	for i in range(3):
		for j in range(i+1, 3):
			o1.plotRhoOffDiag(i, j, 5, lc=colors[i+j-1])
			o2.plotRhoOffDiag(i, j, 5, lc=colors[i+j+2])
	finalizePlot(
		"plots/3nu/rho_offdiag_20_%s.pdf"%fn,
		lloc="lower left",
		)

	for i in range(3):
		for j in range(i+1, 3):
			o1.plotdRhoOffDiag(i, j, 5, lc=colors[i+j-1])
			o2.plotdRhoOffDiag(i, j, 5, lc=colors[i+j+2])
	finalizePlot(
		"plots/3nu/rho_doffdiag_20_%s.pdf"%fn,
		lloc="lower right",
		)

for fn, cl, ol in [
		["ord", (3, 0, 1), (fnoosc20, fNO, fIO)],
		["no", (0, 1), (fNO, fNOd,)],
		["io", (0, 1), (fIO, fIOd,)],
		]:
	for i in range(3):
		for o, c in zip(ol, cl):
			o.plotRhoFin(i, lc=colors[c], ls=styles[i])
	finalizePlot(
		"plots/3nu/rhofin_diag_20_%s.pdf"%fn,
		lloc="upper left",
		xlim=[0, 12],
		ylim=[-0.01, 0.07],
		)
