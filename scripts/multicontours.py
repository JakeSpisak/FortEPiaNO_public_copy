#!/usr/bin/python
import numpy as np
import matplotlib
matplotlib.use("agg")
import matplotlib.pyplot as plt
from dogrid3p1 import Namespace, call_read, fillGrid, indexes, safegetattr

levels=[3., 3.1, 3.3, 3.9]
labels={"e": r"$\alpha=e$", "m": r"$\alpha=\mu$", "t": r"$\alpha=\tau$"}
lstyles=['dashed', 'dashdot', 'dotted']
cmap=matplotlib.cm.get_cmap('CMRmap')
colors=[cmap((l-2.5)/2) for l in levels]
bboxprop = dict(facecolor='white', alpha=0.8, boxstyle='round')


def convert_grid(convert_x, convert_y, pts):
	outpts = np.column_stack((pts, np.zeros(np.shape(pts)[0])))
	outpts = np.column_stack((outpts, np.zeros(np.shape(pts)[0])))
	return outpts


for order in ["no", "io"]:
	fig = plt.figure(figsize=(6,4))
	for ix, grid in enumerate(["e", "m", "t"]):
		args = Namespace(
			gridname="U%s4_%s"%(grid, order),
			verbose=False,
			Neff_ref=0.,
			par_x="U%s4sq"%grid,
			par_y="dm41",
			fix_dm41=False,
			fix_Ue4sq=False if grid=="e" else 0,
			fix_Um4sq=False if grid=="m" else 0,
			fix_Ut4sq=False if grid=="t" else 0,
			)
		mixings, fullgrid, fullobjects = call_read(args)
		fullpoints = list(map(lambda x: safegetattr(x, "Neff", 0.)-args.Neff_ref, fullobjects))
		fullpoints = np.asarray(fullpoints)
		cgrid = {}
		fullgrid = convert_grid(args.par_x, args.par_y, fullgrid)
		for a in indexes.keys():
			if mixings["%s_N"%a] > 1 and getattr(args, "fix_%s"%a) is False:
				cgrid["%s_N"%a] = mixings["%s_N"%a]
				cgrid["%s_min"%a] = mixings["%s_min"%a]
				cgrid["%s_max"%a] = mixings["%s_max"%a]
				cgrid["%s_pts"%a] = mixings["%s_pts"%a]
			else:
				print("fixing %s to %s"%(a, mixings["%s_pts"%a][int(getattr(args, "fix_%s"%a))]))
				cgrid["%s_N"%a] = 1
				cgrid["%s_min"%a] = mixings["%s_pts"%a][int(getattr(args, "fix_%s"%a))]
				cgrid["%s_max"%a] = mixings["%s_pts"%a][int(getattr(args, "fix_%s"%a))]
				cgrid["%s_pts"%a] = np.asarray([mixings["%s_pts"%a][int(getattr(args, "fix_%s"%a))]])
		smallgrid = fillGrid(Namespace(**cgrid))
		smallgrid = convert_grid(args.par_x, args.par_y, smallgrid)
		smallpoints = []
		for fgv, pt in zip(fullgrid, fullpoints):
			for ngv in smallgrid:
				if list(fgv) == list(ngv):
					smallpoints.append(pt)
		smallpoints = np.asarray(smallpoints)
		xv = mixings["%s_pts"%args.par_x]
		yv = mixings["%s_pts"%args.par_y]
		zv = smallpoints.reshape((len(yv), len(xv)))
		plt.contour(
			xv,
			yv,
			zv,
			levels=levels,
			colors=colors,
			linestyles=lstyles[ix],
			extend="both")
		plt.plot(np.nan, label=labels[grid], ls=lstyles[ix], c="k")
	plt.text(1.2e-6, 30, r"$N_{\rm eff}\leq 3.1$", color=colors[1], bbox=bboxprop)
	plt.text(3e-5, 30, r"$N_{\rm eff}\simeq 3.3$", color=colors[2], bbox=bboxprop)
	plt.text(7e-4, 30, r"$N_{\rm eff}\geq 3.9$", color=colors[3], bbox=bboxprop)
	#electron disappearance
	for constr_file in [
			"/home/gariazzo/data/lsn/globalfit/1801.06469/fig4b-contours/cnt-9973-%d.dat"%i
			for i in [1,2,3]]:
		data = np.loadtxt(constr_file)
		plt.plot((1.-np.sqrt(1.-data[:,0]))*0.5, data[:,1], c="r", ls="--")
	#muon disappearance
	if order == "no":
		for constr_file in [
				"/home/gariazzo/data/lsn/minos+/1710.06488v2/minosp_90cl_1710_06488v2_f3.dat"]:
			data = np.loadtxt(constr_file)
			plt.plot(data[:,0], data[:,1], c="b", ls="-.")
	ax=plt.gca()
	ax.tick_params("both", which="both", direction="in",
		left=True, right=True, top=True, bottom=True)
	ax.tick_params("both", which="major", direction="in", length=8)
	plt.xscale("log")
	plt.yscale("log")
	plt.xlabel(r"$|U_{\alpha4}|^2$")
	plt.ylabel(r"$\Delta m^2_{41}$ [eV$^2$]")
	plt.xlim([1e-6, 3e-1])
	plt.ylim([5e-3, 1e2] if order == "no" else [1e-4, 1e2])
	plt.tight_layout(rect=(-0.03, -0.03, 1.02, 1.02))
	plt.legend(loc="lower left")
	plt.savefig("../plots/3p1/angles_%s.pdf"%order)
	plt.close()
