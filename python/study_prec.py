import matplotlib
matplotlib.use("agg")
import matplotlib.pyplot as plt
import numpy as np
from nuDensOutput import colors, markers, styles, finalizePlot, stripRepeated, NuDensRun

labels = {
	"lp": "tol=1e-3",
	"ip": "tol=1e-4",
	"hp": "tol=1e-5",
	"mp": "tol=1e-6",
	}
attr = "Neff"

def plotAttr(attr, ny, ip, lab, folder, color="k", mar=".", nnu=3, prange=(0, 1e3)):
	try:
		obj = NuDensRun(folder, label=lab, nnu=nnu)
	except (IOError, IndexError):
		print("no %s"%lab)
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
			plt.plot(ny + (ip-1)*0.5 - 0.25, ptv, linestyle='none', color=color, marker=mar)

print("3nu")
fig = plt.Figure(figsize=(6,4))
yran = (3.041, 3.046)
for ip, prec in enumerate(["lp", "ip", "hp", "mp"]):
	plt.plot(0, 0, color=colors[ip], marker=".", linestyle='none', label=labels[prec])
	for ny in [20, 40, 70, 100]:
		plotAttr(attr, ny, ip, "Ny=%s %s"%(ny, prec), "OUT/prec_3nu/%s_%s/"%(ny, prec), prange=yran, color=colors[ip])
finalizePlot(
	"plots/prec_3nu.pdf",
	xlab="$N_y$",
	ylab=r"$N_{\rm eff}$",
	xlim=(16, 104),
	ylim=yran,
	# title="3 neutrinos, normal ordering, VLC 2018 best-fit",
	)

print("\n\n3+1")
fig = plt.Figure(figsize=(6,4))
yran = (4.02, 4.08)
for ip, prec in enumerate(["lp", "ip", "hp", "mp"]):
	plt.plot(0, 0, color=colors[ip], marker=".", linestyle='none', label=labels[prec])
	for ny in [15, 20, 25, 30, 40, 70]:
		plotAttr(attr, ny, ip, "Ny=%s %s"%(ny, prec), "OUT/prec_3p1/%s_%s/"%(ny, prec), nnu=4, prange=yran, color=colors[ip])
finalizePlot(
	"plots/prec_3p1.pdf",
	xlab="$N_y$",
	ylab=r"$N_{\rm eff}$",
	xlim=(12, 73),
	ylim=yran,
	# title="3 (NO, VLC 2018) + 1 (Delta=1.29, Ue4=0.01, Um4=Ut4=0) neutrinos",
	lloc="lower right",
	)

print("\n\n3+1 xy")
fig = plt.Figure(figsize=(6,4))
yran = (4.02, 4.08)
for ix, xin in enumerate([1e-3, 9e-4, 8e-4, 5e-4]):
	plt.plot(0, 0, color='k', marker=markers[ix], linestyle='none', label=r"$x_{\rm in}$=%s"%xin)
for iy, nyl in enumerate([2, 3, 4]):
	plt.plot(0, 0, color=colors[iy], marker="o", linestyle='none', label=r"$N_y^{\rm log}$=%s"%nyl)
for iy, nyl in enumerate([2, 3, 4]):
	for ix, xin in enumerate([1e-3, 9e-4, 8e-4, 5e-4]):
		rlab = "Nylog=%s, xin=%s"%(nyl, xin)
		for ny in [20, 25, 30, 35]:
			plotAttr(attr, ny, ix, "Ny=%s %s"%(ny, rlab), "OUT/prec_3p1/%s_hp_%s_%s/"%(ny, nyl, xin),
				nnu=4, prange=yran, color=colors[iy], mar=markers[ix])
finalizePlot(
	"plots/prec_3p1_xy.pdf",
	xlab="$N_y$",
	ylab=r"$N_{\rm eff}$",
	xlim=(17, 38),
	ylim=yran,
	# title="3 (NO, VLC 2018) + 1 (Delta=1.29, Ue4=0.01, Um4=Ut4=0) neutrinos",
	lloc="lower right", legcol=2,
	)
