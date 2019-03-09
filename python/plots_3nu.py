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
fNO40 = NuDensRun("OUT/3nu/40_ci_no/", label="NO")
fIO40 = NuDensRun("OUT/3nu/40_ci_io/", label="IO")
fNO40d = NuDensRun("OUT/3nu/40_df_no/", label="NO damp")
fIO40d = NuDensRun("OUT/3nu/40_df_io/", label="IO damp")
fNO40nm = NuDensRun("OUT/3nu/40_df_nm/", label="NO no mix")
fNO40nd = NuDensRun("OUT/3nu/40_nd_no/", label="NO no damp")
fnoosc40 = NuDensRun("OUT/3nu/40_noosc/", label="no osc Ny=40", full=False)
fNO70 = NuDensRun("OUT/3nu/70_ci_no/", label="NO")
fIO70 = NuDensRun("OUT/3nu/70_ci_io/", label="IO")
fNO70d = NuDensRun("OUT/3nu/70_df_no/", label="NO damp")
fIO70d = NuDensRun("OUT/3nu/70_df_io/", label="IO damp")
fnoosc70 = NuDensRun("OUT/3nu/70_noosc/", label="no osc Ny=70", full=False)
fNO100 = NuDensRun("OUT/3nu/100_ci_no/", label="NO")
fIO100 = NuDensRun("OUT/3nu/100_ci_io/", label="IO")
fNO100d = NuDensRun("OUT/3nu/100_df_no/", label="NO damp")
fIO100d = NuDensRun("OUT/3nu/100_df_io/", label="IO damp")
fnoosc100 = NuDensRun("OUT/3nu/100_noosc/", label="no osc Ny=100", full=False)

for i in range(3):
	fnoosc20.plotRhoDiag(i, 4, styles[i], lc=colors[0])
	fnoosc40.plotRhoDiag(i, 9, styles[i], lc=colors[1])
	fnoosc70.plotRhoDiag(i, 17, styles[i], lc=colors[2])
	fnoosc100.plotRhoDiag(i, 20, styles[i], lc=colors[4])
finalizePlot(
	"plots/3nu/rho_diag_noosc.pdf",
	xlab="$x$",
	ylab=r"$\rho/\rho_{eq}-1$",
	xscale="log",
	lloc="upper left",
	)

for i in range(3):
	fnoosc20.plotRhoDiagY(i, 5., styles[i], lc=colors[0])
	fnoosc40.plotRhoDiagY(i, 5., styles[i], lc=colors[1])
	fnoosc70.plotRhoDiagY(i, 5., styles[i], lc=colors[2])
	fnoosc100.plotRhoDiagY(i, 5., styles[i], lc=colors[4])
finalizePlot(
	"plots/3nu/rho_diag_noosc_y.pdf",
	title="$y=5$",
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

fNOd.label = "NO Ny=20"
fIOd.label = "IO Ny=20"
fNO100d.label = "NO Ny=100"
fIO100d.label = "IO Ny=100"
for i in range(3):
	fnoosc100.plotRhoDiagY(i, 5, styles[i], lc=colors[3])
	fNOd.plotRhoDiagY(i, 5, styles[i], lc=colors[0])
	fIOd.plotRhoDiagY(i, 5, styles[i], lc=colors[1])
	fNO100d.plotRhoDiagY(i, 5, styles[i], lc=colors[2])
	fIO100d.plotRhoDiagY(i, 5, styles[i], lc=colors[4])
finalizePlot(
	"plots/3nu/rho_diag_Ny_damp.pdf",
	title="with damping terms only",
	xlab="$x$",
	ylab=r"$\rho/\rho_{eq}-1$",
	xscale="log",
	lloc="upper left",
	)

for fn, ol in [
		["ord", [fNOd, fIOd, fNO100d, fIO100d]],
		]:
	for i in range(3):
		for j in range(i+1, 3):
			for io, o in enumerate(ol):
				o.plotRhoOffDiagY(i, j, 5, lc=colors[io], ls=styles[i+j], im=False)
	finalizePlot(
		"plots/3nu/rho_offdiag_Ny_%s.pdf"%fn,
		lloc="lower left",
		)

	for i in range(3):
		for j in range(i+1, 3):
			for io, o in enumerate(ol):
				o.plotdRhoOffDiagY(i, j, 5, lc=colors[io], ls=styles[i+j], im=False)
	finalizePlot(
		"plots/3nu/rho_doffdiag_Ny_%s.pdf"%fn,
		lloc="lower right",
		)

fNOd.plotDeltaZ(fnoosc100, lc=colors[0])
fIOd.plotDeltaZ(fnoosc100, lc=colors[1])
fNO100d.plotDeltaZ(fnoosc100, lc=colors[2])
fIO100d.plotDeltaZ(fnoosc100, lc=colors[4])
finalizePlot(
	"plots/3nu/deltaz_damp.pdf",
	title=r"$\Delta z$ with respect to Ny=100, no oscillations",
	xlab="$x$",
	ylab=r"$z-z_{\rm no osc}$",
	xscale="log",
	)
