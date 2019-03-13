import sys
import os
import argparse
import ast
import glob
import re
import numpy as np
from scipy.interpolate import interpn
import prepareIni
import matplotlib
matplotlib.use("agg")
import matplotlib.pyplot as plt
from nuDensOutput import colors, markers, styles, finalizePlot, stripRepeated, NuDensRun


cmap=matplotlib.cm.get_cmap('CMRmap')
labels = {
	"dm41": r"$\Delta m^2_{41}$ [eV$^2$]",
	"ssq14": r"$\sin^2\theta_{14}$",
	"ssq24": r"$\sin^2\theta_{24}$",
	"ssq34": r"$\sin^2\theta_{34}$",
	}
indexes = {"dm41": 0, "ssq14": 1, "ssq24": 2, "ssq34": 3}


class Namespace:
	def __init__(self, **kwargs):
		self.__dict__.update(kwargs)


def safegetattr(obj, attr, default):
	try:
		value = getattr(obj, attr)
	except AttributeError:
		return default
	else:
		return value


def setParser():
	parser = argparse.ArgumentParser(prog='dogrid3p1.py')
	parser.add_argument(
		'gridname',
		help='the name of the grid you want to use'
		)
	subparsers = parser.add_subparsers(
		help='sub-command help',
		dest='cmd',
		)
	subparsers.required = True

	parser_plot = subparsers.add_parser(
		'plot',
		help='plot some output'
		)
	parser_plot.add_argument(
		'--Neff_active',
		type=float,
		default=3.044,
		help='reference value of Neff when only active neutrinos are considered',
		)
	parser_plot.add_argument(
		'--par_x',
		choices=["ssq14", "ssq24", "ssq34"],
		default="ssq24",
		help='parameter for the x axis',
		)
	parser_plot.add_argument(
		'--par_y',
		choices=["dm41", "ssq14", "ssq24", "ssq34"],
		default="dm41",
		help='parameter for the y axis',
		)
	parser_plot.add_argument(
		'--fix_dm41',
		default=False,
		help='when requesting a subplot, fix dm41 to the given index (not value!)',
		)
	parser_plot.add_argument(
		'--fix_ssq14',
		default=0,
		help='when requesting a subplot, fix ssq14 to the given index (not value!)',
		)
	parser_plot.add_argument(
		'--fix_ssq24',
		default=False,
		help='when requesting a subplot, fix ssq24 to the given index (not value!)',
		)
	parser_plot.add_argument(
		'--fix_ssq34',
		default=0,
		help='when requesting a subplot, fix ssq34 to the given index (not value!)',
		)
	parser_plot.add_argument(
		'--title',
		default="",
		help='title for the plot',
		)
	parser_plot.set_defaults(func=call_plot)

	parser_prepare = subparsers.add_parser(
		'set',
		help='set the grid and create ini files'
		)
	parser_prepare.add_argument(
		'--dm41_min',
		type=float,
		default=1e-2,
		help='minimum value of Delta m^2_{41}',
		)
	parser_prepare.add_argument(
		'--dm41_max',
		type=float,
		default=1e2,
		help='maximum value of Delta m^2_{41}',
		)
	parser_prepare.add_argument(
		'--dm41_N',
		type=int,
		default=5,
		help='number of points in Delta m^2_{41}',
		)
	parser_prepare.add_argument(
		'--ssq14_min',
		type=float,
		default=1e-4,
		help='minimum value of sin^2 theta_{14}',
		)
	parser_prepare.add_argument(
		'--ssq14_max',
		type=float,
		default=1,
		help='maximum value of sin^2 theta_{14}',
		)
	parser_prepare.add_argument(
		'--ssq14_N',
		type=int,
		default=5,
		help='number of points in sin^2 theta_{14}',
		)
	parser_prepare.add_argument(
		'--ssq24_min',
		type=float,
		default=1e-4,
		help='minimum value of sin^2 theta_{24}',
		)
	parser_prepare.add_argument(
		'--ssq24_max',
		type=float,
		default=1,
		help='maximum value of sin^2 theta_{24}',
		)
	parser_prepare.add_argument(
		'--ssq24_N',
		type=int,
		default=5,
		help='number of points in sin^2 theta_{24}',
		)
	parser_prepare.add_argument(
		'--ssq34_min',
		type=float,
		default=1e-4,
		help='minimum value of sin^2 theta_{34}',
		)
	parser_prepare.add_argument(
		'--ssq34_max',
		type=float,
		default=1,
		help='maximum value of sin^2 theta_{34}',
		)
	parser_prepare.add_argument(
		'--ssq34_N',
		type=int,
		default=5,
		help='number of points in sin^2 theta_{34}',
		)
	parser_prepare.set_defaults(func=call_prepare)

	parser_read = subparsers.add_parser(
		'read',
		help='read the output files and print resume'
		)
	parser_read.set_defaults(func=call_read)

	parser_run = subparsers.add_parser(
		'run',
		help='submit the jobs in the grid'
		)
	parser_run.set_defaults(func=call_run)
	return parser


def fillGrid(args):
	return 10**(np.mgrid[
		np.log10(args.dm41_min):np.log10(args.dm41_max):args.dm41_N*1j,
		np.log10(args.ssq14_min):np.log10(args.ssq14_max):args.ssq14_N*1j,
		np.log10(args.ssq24_min):np.log10(args.ssq24_max):args.ssq24_N*1j,
		np.log10(args.ssq34_min):np.log10(args.ssq34_max):args.ssq34_N*1j,
		]).reshape(4, args.dm41_N*args.ssq14_N*args.ssq24_N*args.ssq34_N).T


def write_grid_cfg(args):
	with open("grids/%s/params.cfg"%args.gridname, "w") as _f:
		for a in ["dm41", "ssq14", "ssq24", "ssq34"]:
			if "ssq" in a:
				_f.write("%s: min=%s, max=%s, N=%s, ssq2th=%s\n"%(
					a,
					getattr(args, "%s_min"%a),
					getattr(args, "%s_max"%a),
					getattr(args, "%s_N"%a),
					getattr(args, a.replace("ssq", "sinsq2th")),
					))
			else:
				_f.write("%s: min=%s, max=%s, N=%s\n"%(
					a,
					getattr(args, "%s_min"%a),
					getattr(args, "%s_max"%a),
					getattr(args, "%s_N"%a),
					))


def read_grid_cfg(gridname):
	with open("grids/%s/params.cfg"%args.gridname) as _f:
		text = _f.readlines()
	values = {}
	for a in ["dm41", "ssq14", "ssq24", "ssq34"]:
		for l in text:
			res = re.match("%s: min=([0-9\.e\+\-]+), max=([0-9\.e\+\-]+), N=([0-9]+)\n"%a, l)
			if res:
				try:
					values["%s_min"%a] = float(res.group(1))
					values["%s_max"%a] = float(res.group(2))
					values["%s_N"%a] = int(res.group(3))
				except (IndexError, TypeError, ValueError):
					print("cannot read line! %s"%l)
				else:
					if values["%s_min"%a] == 0 or values["%s_max"%a] == 0:
						values["%s_pts"%a] = np.linspace(
							values["%s_min"%a],
							values["%s_max"%a],
							values["%s_N"%a])
					else:
						values["%s_pts"%a] = np.logspace(
							np.log10(values["%s_min"%a]),
							np.log10(values["%s_max"%a]),
							values["%s_N"%a])
			else:
				res = re.match("%s: min=([0-9\.e\+\-]+), max=([0-9\.e\+\-]+), N=([0-9]+), ssq2th=(True|False)\n"%a, l)
				if res:
					try:
						values["%s_min"%a] = float(res.group(1))
						values["%s_max"%a] = float(res.group(2))
						values["%s_N"%a] = int(res.group(3))
						values[a.replace("ssq", "sinsq2th")] = ast.literal_eval(res.group(4))
					except (IndexError, TypeError, ValueError):
						print("cannot read line! %s"%l)
					else:
						if values["%s_min"%a] == 0 or values["%s_max"%a] == 0:
							values["%s_pts"%a] = np.linspace(
								values["%s_min"%a],
								values["%s_max"%a],
								values["%s_N"%a])
						else:
							values["%s_pts"%a] = np.logspace(
								np.log10(values["%s_min"%a]),
								np.log10(values["%s_max"%a]),
								values["%s_N"%a])
	for a in ["ssq14", "ssq24", "ssq34"]:
		if a.replace("ssq", "sinsq2th") not in values.keys():
			values[a.replace("ssq", "sinsq2th")] = False
	return values, fillGrid(Namespace(**values))


def contourplot(xv, yv, values, points,
		fname=None,
		levels=[0., 0.1, 0.3, 0.5, 0.7, 0.9],
		title=None,
		xlab=None, ylab=None,
		xlim=None, ylim=None,
		):
	zv = points.reshape((len(xv), len(yv)))
	fig = plt.Figure(figsize=(6,4))
	cf = plt.contourf(xv, yv, zv, levels=levels, cmap=matplotlib.cm.get_cmap('CMRmap'), extend="max")
	cf.cmap.set_over(cmap(0.95))
	cbar = plt.colorbar(cf)
	cbar.ax.set_ylabel(r"$\Delta N_{\rm eff}$")
	if title is not None:
		plt.title(title)
	ax=plt.gca()
	ax.tick_params("both", which="both", direction="out",
		left=True, right=True, top=True, bottom=True)
	ax.tick_params("both", which="major", direction="out", length=8)
	plt.xscale("log")
	plt.yscale("log")
	if xlab is not None:
		plt.xlabel(xlab)
	if ylab is not None:
		plt.ylabel(ylab)
	if xlim is not None:
		plt.xlim(xlim)
	if ylim is not None:
		plt.ylim(ylim)
	plt.tight_layout()
	if fname is not None:
		plt.savefig(fname)
	plt.close()


def call_plot(args):
	mixings, fullgrid, fullobjects = call_read(args)
	fullpoints = list(map(lambda x: safegetattr(x, "Neff", 0.)-args.Neff_active, fullobjects))
	fullpoints = np.asarray(fullpoints)
	cgrid = {}
	for a in indexes.keys():
		if args.par_x == a or args.par_y == a:
			setattr(args, "fix_%s"%a, False)
		if mixings["%s_N"%a] > 1 and getattr(args, "fix_%s"%a) is False:
			cgrid["%s_N"%a] = mixings["%s_N"%a]
			cgrid["%s_min"%a] = mixings["%s_min"%a]
			cgrid["%s_max"%a] = mixings["%s_max"%a]
			cgrid["%s_pts"%a] = mixings["%s_pts"%a]
		else:
			cgrid["%s_N"%a] = 1
			cgrid["%s_min"%a] = mixings["%s_pts"%a][int(getattr(args, "fix_%s"%a))]
			cgrid["%s_max"%a] = mixings["%s_pts"%a][int(getattr(args, "fix_%s"%a))]
			cgrid["%s_pts"%a] = np.asarray([mixings["%s_pts"%a][int(getattr(args, "fix_%s"%a))]])
	smallgrid = fillGrid(Namespace(**cgrid))
	xv = mixings["%s_pts"%args.par_x]
	yv = mixings["%s_pts"%args.par_y]
	smallpoints = []
	for [p1, p2, p3, p4], pt in zip(fullgrid, fullpoints):
		for n1, n2, n3, n4 in smallgrid:
			if (p1==n1 and p2==n2 and p3==n3 and p4==n4):
				smallpoints.append(pt)
	smallpoints = np.asarray(smallpoints)
	if "ssq" in args.par_x and mixings[args.par_x.replace("ssq", "sinsq2th")]:
		xlab = labels[args.par_x].replace(r"\theta", r"2\theta")
	else:
		xlab = labels[args.par_x]
	if "ssq" in args.par_y and mixings[args.par_y.replace("ssq", "sinsq2th")]:
		ylab = labels[args.par_y].replace(r"\theta", r"2\theta")
	else:
		ylab = labels[args.par_y]
	contourplot(
		xv, yv, mixings, smallpoints,
		fname="grids/%s/plots/%s_%s.pdf"%(args.gridname, args.par_x, args.par_y),
		xlab=xlab, ylab=ylab, title=args.title,
		)


def call_prepare(args):
	if not os.path.exists("grids/%s/ini/"%args.gridname):
		os.makedirs("grids/%s/ini/"%args.gridname)
	if not os.path.exists("grids/%s/OUT/"%args.gridname):
		os.makedirs("grids/%s/OUT/"%args.gridname)
	if not os.path.exists("grids/%s/plots/"%args.gridname):
		os.makedirs("grids/%s/plots/"%args.gridname)
	files = list(glob.iglob("grids/%s/ini/*.ini"%args.gridname))
	list(map(lambda x: os.remove(x), files))
	for a in ["dm41", "ssq14", "ssq24", "ssq34"]:
		if getattr(args, "%s_N"%a) == 1:
			pmin = getattr(args, "%s_min"%a)
			pmax = getattr(args, "%s_max"%a)
			if pmin == 0. or pmax == 0.:
				setattr(args, "%s_min"%a,
					0.5 * (pmin + pmax)
					)
			else:
				setattr(args, "%s_min"%a,
					10**(0.5 *
						(np.log10(pmin) + np.log10(pmax))
					))
			setattr(args, "%s_max"%a, getattr(args, "%s_min"%a))
	write_grid_cfg(args)
	grid = fillGrid(args)
	for dm41, ssq14, ssq24, ssq34 in grid:
		if args.sinsq2th14:
			ssq14 = 0.5*(1.-np.sqrt(1.-ssq14))
		if args.sinsq2th24:
			ssq24 = 0.5*(1.-np.sqrt(1.-ssq24))
		if args.sinsq2th34:
			ssq34 = 0.5*(1.-np.sqrt(1.-ssq34))
		prep = [
			"grids/%s/ini/%s_%s_%s_%s.ini"%(args.gridname, dm41, ssq14, ssq24, ssq34),
			"grids/%s/OUT/%s_%s_%s_%s/"%(args.gridname, dm41, ssq14, ssq24, ssq34),
			"3+1",
			"damping",
			"--dlsoda_rtol=1e-5",
			"--dlsoda_atol=1e-5",
			"--Nx=2000",
			"--Ny=30",
			"--Nylog=4",
			"--y_cen=1",
			"--x_in=0.0005",
			"--default_sterile=None",
			"--dm41=%s"%dm41,
			"--th14=%s"%ssq14,
			"--th24=%s"%ssq24,
			"--th34=%s"%ssq34,
			]
		parser = prepareIni.setParser()
		rargs = parser.parse_args(prep)
		values = prepareIni.getIniValues(rargs)
		prepareIni.writeIni(rargs.inifile, values)


def call_read(args):
	values, grid = read_grid_cfg(args.gridname)
	objects = []
	for dm41, ssq14, ssq24, ssq34 in grid:
		lab = (r"dm41=%s "%dm41
			+ r"ssq14=%s "%ssq14
			+ r"ssq24=%s "%ssq24
			+ r"ssq34=%s "%ssq34
			)
		if values["sinsq2th14"]:
			ssq14 = 0.5*(1.-np.sqrt(1.-ssq14))
		if values["sinsq2th24"]:
			ssq24 = 0.5*(1.-np.sqrt(1.-ssq24))
		if values["sinsq2th34"]:
			ssq34 = 0.5*(1.-np.sqrt(1.-ssq34))
		folder = "grids/%s/OUT/%s_%s_%s_%s/"%(args.gridname, dm41, ssq14, ssq24, ssq34)
		obj = None
		try:
			obj = NuDensRun(folder, label=lab, nnu=4, rho=False)
		except (IOError, IndexError):
			print("no %s"%lab)
		objects.append(obj)
	return values, grid, objects


def call_run(args):
	files = list(glob.iglob("grids/%s/ini/*.ini"%args.gridname))
	print("submitting the grid %s"%args.gridname)
	for f in files:
		os.system(
			"clusterlauncher -N r%s -n 1 --openmp -q short-seq bin/nuDens.exe %s"%(
				f.split(os.sep)[-1].replace(".ini", ""),
				f,
				))


if __name__=='__main__':
	parser = setParser()
	args = parser.parse_args(sys.argv[1:])
	print(args)
	args.func(args)
