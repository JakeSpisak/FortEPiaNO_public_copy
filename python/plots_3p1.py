import traceback
import matplotlib
matplotlib.use("agg")
import matplotlib.pyplot as plt
import numpy as np
from nuDensOutput import colors, styles, finalizePlot, stripRepeated, NuDensRun

fNO40d = NuDensRun("OUT/3nu/40_df_no/", label="3nu NO damp")
s3p1_df_lp = NuDensRun("OUT/prec_3p1/15_lp/", nnu=4, label="3p1 low prec, Ny=15")
s3p1_df_hp = NuDensRun("OUT/prec_3p1/20_mp/", nnu=4, label="3p1 high prec, Ny=20")
s3p1_df_hp40 = NuDensRun("OUT/prec_3p1/40_lp/", nnu=4, label="3p1 low prec, Ny=40")
s3p1_df_hp40 = NuDensRun("OUT/prec_3p1/40_mp/", nnu=4, label="3p1 high prec, Ny=40")
s3p1_df_hp70 = NuDensRun("OUT/prec_3p1/70_lp/", nnu=4, label="3p1 low prec, Ny=70")
s3p1_df_hp70 = NuDensRun("OUT/prec_3p1/70_mp/", nnu=4, label="3p1 high prec, Ny=70", plots=True)

for i in range(3):
	fNO40d.plotRhoDiagY(i, 5, styles[i], lc=colors[3])
for i in range(4):
	s3p1_df_lp.plotRhoDiagY(i, 5, styles[i], lc=colors[0])
for i in range(4):
	s3p1_df_hp.plotRhoDiagY(i, 5, styles[i], lc=colors[1])
finalizePlot(
	"plots/3p1/rho_diag_prec.pdf",
	xlab="$x$",
	ylab=r"$\rho/\rho_{eq}-1$",
	xscale="log",
	lloc="lower right",
	)

fNO40d.plotZ(lc=colors[3])
s3p1_df_lp.plotZ(lc=colors[0])
s3p1_df_hp.plotZ(lc=colors[1])
finalizePlot(
	"plots/3p1/z_prec.pdf",
	xlab="$x$",
	ylab=r"$z$",
	xscale="log",
	)

# s3p1_df_lp.plotDeltaZ(fNO40d, lc=colors[0])
# s3p1_df_hp.plotDeltaZ(fNO40d, lc=colors[1])
# finalizePlot(
	# "plots/3p1/deltaz_prec.pdf",
	# xlab="$x$",
	# ylab=r"$z-z_{\rm 3nu}$",
	# xscale="log",
	# )

for ol, fn in [
		[[fNO40d, s3p1_df_lp, s3p1_df_hp], "prec"],
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
		["ord", (3, 0, 1), (fNO40d, s3p1_df_lp, s3p1_df_hp)],
		]:
	for o, c in zip(ol, cl):
		for i in range(4):
			try:
				o.plotRhoFin(i, ls=styles[i], lc=colors[c])
			except IndexError:
				pass#print(traceback.format_exc())
	finalizePlot(
		"plots/3p1/rhofin_diag_%s.pdf"%fn,
		lloc="lower left",
		# xlim=[0, 12],
		# ylim=[-0.01, 0.07],
		)
