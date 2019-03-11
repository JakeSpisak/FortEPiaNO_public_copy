import matplotlib
matplotlib.use("agg")
import matplotlib.pyplot as plt
import numpy as np
from nuDensOutput import colors, styles, finalizePlot, stripRepeated, NuDensRun

labels = {
	"lp": "tol=1e-3",
	"ip": "tol=1e-4",
	"hp": "tol=1e-5",
	"mp": "tol=1e-6",
	}
attr = "Neff"

def plotAttr(attr, ny, ip, prec, folder, nnu=3, prange=(0, 1e3)):
	try:
		obj = NuDensRun("%s/%s_%s/"%(folder, ny, prec), label="Ny=%s %s"%(ny, prec), nnu=nnu)
	except (IOError, IndexError):
		print("no Ny=%s %s"%(ny, prec))
	else:
		try:
			ptv = getattr(obj, attr)
		except AttributeError:
			pass
		else:
			if ptv < prange[0]:
				mar = r'$\downarrow$'
				ptv = prange[0] + (prange[1] - prange[0]) * 0.02
			elif ptv > prange[1]:
				mar = r'$\uparrow$'
				ptv = prange[1] - (prange[1] - prange[0]) * 0.02
			else:
				mar = "."
			plt.plot(ny + (ip-1)*0.5 - 0.25, ptv, linestyle='none', color="%s"%colors[ip], marker=mar)

print("3nu")
fig = plt.Figure(figsize=(6,4))
addlabel = True
yran = (3.041, 3.046)
for ip, prec in enumerate(["lp", "ip", "hp", "mp"]):
	plt.plot(0, 0, "%s."%colors[ip], label=labels[prec])
	for ny in [20, 40, 70, 100]:
		plotAttr(attr, ny, ip, prec, "OUT/prec_3nu/", prange=yran)
finalizePlot(
	"plots/prec_3nu.pdf",
	xlab="$N_y$",
	ylab=r"$N_{\rm eff}$",
	xlim=(16, 104),
	ylim=yran,
	title="3 neutrinos, normal ordering, VLC 2018 best-fit",
	)

print("\n\n3+1")
fig = plt.Figure(figsize=(6,4))
addlabel = True
yran = (4.02, 4.08)
for ip, prec in enumerate(["lp", "ip", "hp", "mp"]):
	plt.plot(0, 0, "%s."%colors[ip], label=labels[prec])
	for ny in [15, 20, 25, 30, 40, 70]:
		plotAttr(attr, ny, ip, prec, "OUT/prec_3p1/", nnu=4, prange=yran)
	addlabel = False
finalizePlot(
	"plots/prec_3p1.pdf",
	xlab="$N_y$",
	ylab=r"$N_{\rm eff}$",
	xlim=(12, 73),
	ylim=yran,
	title="3 (NO, VLC 2018) + 1 (Delta=1.29, Ue4=0.01, Um4=Ut4=0) neutrinos",
	lloc="lower right",
	)
