import matplotlib
matplotlib.use("agg")
import matplotlib.pyplot as plt
import numpy as np
from nuDensOutput import colors, styles, stripRepeated, NuDensRun

fNO = NuDensRun("OUT/3nu/20_ci_no/", label="NO")
fIO = NuDensRun("OUT/3nu/20_ci_io/", label="IO")
fNOd = NuDensRun("OUT/3nu/20_df_no/", label="NO damp")
fIOd = NuDensRun("OUT/3nu/20_df_io/", label="IO damp")
fnoosc20 = NuDensRun("OUT/3nu/20_noosc/", label="no osc Ny=20", full=False)
fnoosc40 = NuDensRun("OUT/3nu/40_noosc/", label="no osc Ny=40", full=False)

for i in range(3):
    fnoosc20.plotRhoDiag(i, 5, styles[i], lc=colors[3])
    fNO.plotRhoDiag(i, 5, styles[i], lc=colors[0])
    fIO.plotRhoDiag(i, 5, styles[i], lc=colors[1])
plt.xscale("log")
plt.xlabel("$x$")
plt.ylabel(r"$\rho/\rho_{eq}-1$")
plt.legend(loc="upper left")
plt.tight_layout()
plt.savefig("plots/3nu/rho_diag_20_ord.pdf")
plt.close()

for i in range(3):
	fNO.plotRhoDiag(i, 5, styles[i], lc=colors[0])
	fIO.plotRhoDiag(i, 5, styles[i], lc=colors[1])
	fNOd.plotRhoDiag(i, 5, styles[i], lc=colors[2])
	fIOd.plotRhoDiag(i, 5, styles[i], lc=colors[4])
plt.xscale("log")
plt.xlabel("$x$")
plt.ylabel(r"$\rho/\rho_{eq}-1$")
plt.legend(loc="upper left")
plt.tight_layout()
plt.savefig("plots/3nu/rho_diag_20_damp.pdf")
plt.close()

for i in range(3):
	fnoosc20.plotRhoDiag(i, 3, styles[i], lc=colors[3])
	fnoosc40.plotRhoDiag(i, 8, styles[i], lc=colors[5])
plt.xscale("log")
plt.xlabel("$x$")
plt.ylabel(r"$\rho/\rho_{eq}-1$")
plt.legend(loc="upper left")
plt.tight_layout()
plt.savefig("plots/3nu/rho_diag_noosc.pdf")
plt.close()

fnoosc20.plotZ(lc=colors[3])
fNO.plotZ(lc=colors[0])
fIO.plotZ(lc=colors[1])
fNOd.plotZ(lc=colors[2])
fIOd.plotZ(lc=colors[4])
plt.xscale("log")
plt.xlabel("$x$")
plt.ylabel(r"$z$")
plt.legend()
plt.tight_layout()
plt.savefig("plots/3nu/z_20.pdf")
plt.close()

fNO.plotDeltaZ(fnoosc20, lc=colors[0])
fIO.plotDeltaZ(fnoosc20, lc=colors[1])
fNOd.plotDeltaZ(fnoosc20, lc=colors[2])
fIOd.plotDeltaZ(fnoosc20, lc=colors[4])
plt.xscale("log")
plt.xlabel("$x$")
plt.ylabel(r"$z$")
plt.legend()
plt.tight_layout()
plt.savefig("plots/3nu/deltaz_20.pdf")
plt.close()

for o1, o2, fn in [
		[fNO, fIO, "ord"],
		[fNO, fNOd, "no"],
		[fIO, fIOd, "io"],
		]:
	for i in range(3):
		for j in range(i+1, 3):
			o1.plotRhoOffDiag(i, j, 5, lc=colors[i+j-1])
			o2.plotRhoOffDiag(i, j, 5, lc=colors[i+j+2])
	plt.legend(loc="lower left")
	plt.tight_layout()
	plt.savefig("plots/3nu/rho_offdiag_20_%s.pdf"%fn)
	plt.close()

	for i in range(3):
		for j in range(i+1, 3):
			o1.plotdRhoOffDiag(i, j, 5, lc=colors[i+j-1])
			o2.plotdRhoOffDiag(i, j, 5, lc=colors[i+j+2])
	plt.legend(loc="lower right")
	plt.tight_layout()
	plt.savefig("plots/3nu/rho_doffdiag_20_%s.pdf"%fn)
	plt.close()
