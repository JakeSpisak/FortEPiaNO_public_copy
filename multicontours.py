#!/usr/bin/python
import numpy as np
import matplotlib
matplotlib.use("agg")
import matplotlib.pyplot as plt
from dogrid3p1 import Namespace, call_read, fillGrid, indexes, safegetattr

levels=[3., 3.1, 3.3, 3.5, 3.7, 3.9]
labels={"e": r"$\alpha=e$", "m": r"$\alpha=\mu$", "t": r"$\alpha=\tau$"}
lstyles=['dashed', 'dashdot', 'dotted']
cmap=matplotlib.cm.get_cmap('CMRmap')
colors=[cmap(l-3) for l in levels]


def convert_grid(convert_x, convert_y, pts):
	outpts = np.column_stack((pts, np.zeros(np.shape(pts)[0])))
	outpts = np.column_stack((outpts, np.zeros(np.shape(pts)[0])))
	return outpts


for order in ["no", "io"]:
	fig = plt.figure(figsize=(6,4))
	cf = plt.contourf(
		[np.nan, np.nan],
		[np.nan, np.nan],
		[[np.nan, np.nan], [np.nan, np.nan]],
		levels=levels,
		colors=colors,
		extend="both")
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
	cf.cmap.set_over(cmap(0.95))
	cf.cmap.set_under(cmap(0.0))
	cbar = plt.colorbar(cf, filled=True)
	cbar.ax.set_ylabel(r"$N_{\rm eff}$")
	ax=plt.gca()
	ax.tick_params("both", which="both", direction="in",
		left=True, right=True, top=True, bottom=True)
	ax.tick_params("both", which="major", direction="in", length=8)
	plt.xscale("log")
	plt.yscale("log")
	plt.xlabel(r"$|U_{\alpha4}|^2$")
	plt.ylabel(r"$\Delta m^2_{41}$ [eV$^2$]")
	plt.xlim([1e-6, 3e-1])
	plt.ylim([1e-5, 1e2])
	plt.tight_layout(rect=(-0.03, -0.05, 1.05, 1.02))
	plt.legend()
	plt.savefig("plots/3p1/angles_%s.pdf"%order)
	plt.close()
