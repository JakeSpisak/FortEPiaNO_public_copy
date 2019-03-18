import sys
import os
import argparse
import ast
import glob
import re
import numpy as np
from scipy.interpolate import interpn
import shutil
import subprocess
import matplotlib
matplotlib.use("agg")
import matplotlib.pyplot as plt
import prepareIni
from nuDensOutput import colors, markers, styles, finalizePlot, stripRepeated, NuDensRun


cmap=matplotlib.cm.get_cmap('CMRmap')
labels = {
	"dm41": r"$\Delta m^2_{41}$ [eV$^2$]",
	"Ue4sq": r"$|U_{e4}|^2$",
	"Um4sq": r"$|U_{\mu4}|^2$",
	"Ut4sq": r"$|U_{\tau4}|^2$",
	"sinsqth14": r"$\sin^2\theta_{14}$",
	"sinsqth24": r"$\sin^2\theta_{24}$",
	"sinsqth34": r"$\sin^2\theta_{34}$",
	"sinsq2th_ee": r"$\sin^22\theta_{ee}$",
	"sinsq2th_mumu": r"$\sin^22\theta_{\mu\mu}$",
	"sinsq2th_emu": r"$\sin^22\theta_{e\mu}$",
	}
indexes = {"dm41": 0, "Ue4sq": 1, "Um4sq": 2, "Ut4sq": 3}


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
		default=3.043,
		help='reference value of Neff when only active neutrinos are considered',
		)
	parser_plot.add_argument(
		'--par_x',
		choices=["Ue4sq", "Um4sq", "Ut4sq", "sinsqth14", "sinsqth24", "sinsqth34", "sinsq2th_ee", "sinsq2th_mumu", "sinsq2th_emu"],
		default="Um4sq",
		help='parameter for the x axis',
		)
	parser_plot.add_argument(
		'--par_y',
		choices=["dm41", "Ue4sq", "Um4sq", "Ut4sq", "sinsqth14", "sinsqth24", "sinsqth34", "sinsq2th_ee", "sinsq2th_mumu", "sinsq2th_emu"],
		default="dm41",
		help='parameter for the y axis',
		)
	for a, v in [
			["dm41", False],
			["Ue4sq", 0],
			["Um4sq", False],
			["Ut4sq", 0],
			]:
		parser_plot.add_argument(
			'--fix_%s'%a,
			default=v,
			help='when requesting a subplot, fix %s to the given index (not value!)'%a,
			)
	parser_plot.add_argument(
		'--title',
		default="",
		help='title for the plot',
		)
	parser_plot.add_argument(
		'--filename',
		default="",
		help='name of the file where to save the plot',
		)
	parser_plot.add_argument(
		'--bestfit_x',
		type=float,
		default=-1,
		help='position of the best-fit in x. Set to negative to ignore',
		)
	parser_plot.add_argument(
		'--bestfit_y',
		type=float,
		default=-1,
		help='position of the best-fit in y. Set to negative to ignore',
		)
	parser_plot.add_argument(
		'--bestfit_upper',
		action="store_true",
		help='The angle specified in the bestfit is only an upper limit',
		)
	parser_plot.add_argument(
		'--lsn_contours',
		action="store_true",
		help='Plot contours from global fit',
		)
	parser_plot.set_defaults(func=call_plot)

	parser_prepare = subparsers.add_parser(
		'set',
		help='set the grid and create ini files'
		)
	parser_prepare.add_argument(
		'--Nx',
		type=int,
		default=2000,
		help='number of points in x where to save the output',
		)
	parser_prepare.add_argument(
		'--Ny',
		type=int,
		default=30,
		help='number of momenta',
		)
	parser_prepare.add_argument(
		'--Nylog',
		type=int,
		default=4,
		help='number of log-spaced points in 0.01 < y < 1',
		)
	parser_prepare.add_argument(
		'--x_in',
		type=float,
		default=5e-4,
		help='initial value in x',
		)
	parser_prepare.add_argument(
		'--tolerance',
		type=float,
		default=1e-5,
		help='tolerance for DLSODA',
		)
	for a, l, mi, ma, n in [
			["dm41", r"\Delta m^2_{41}", 1e-2, 1e2, 5],
			["Ue4sq", r"|U_{e4}|^2", 1e-6, 0.1, 6],
			["Um4sq", r"|U_{\mu4}|^2", 1e-6, 0.1, 6],
			["Ut4sq", r"|U_{\tau4}|^2", 1e-6, 0.1, 6],
			]:
		parser_prepare.add_argument(
			'--%s_min'%a,
			type=float,
			default=mi,
			help=r'minimum value of %s'%a,
			)
		parser_prepare.add_argument(
			'--%s_max'%a,
			type=float,
			default=ma,
			help=r'maximum value of %s'%a,
			)
		parser_prepare.add_argument(
			'--%s_N'%a,
			type=int,
			default=n,
			help=r'number of points in %s'%a,
			)
	parser_prepare.add_argument(
		'--ordering',
		choices=["NO", "IO"],
		default="NO",
		help='define the mass ordering for the three active neutrinos'
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
	parser_run.add_argument(
		'-f',
		'--failed_only',
		action="store_true",
		help='resubmit only failed or incomplete runs',
		)
	parser_run.add_argument(
		'-r',
		'--remove_existing',
		action="store_true",
		help='remove the folder before resubmitting the failed run',
		)
	parser_run.add_argument(
		'-s',
		'--first-index',
		type=int,
		default=0,
		help='starting index for the runs to submit (if there are a lot of runs)',
		)
	parser_run.add_argument(
		'-t',
		'--last-index',
		type=int,
		default=1000,
		help='ending index for the runs to submit (if there are a lot of runs)',
		)
	parser_run.add_argument(
		'--expected_lines',
		type=int,
		default=2000,
		help=r'number of expected lines (values in x) in output',
		)
	parser_run.set_defaults(func=call_run)
	return parser


def fillGrid(args):
	return 10**(np.mgrid[
		np.log10(args.dm41_min):np.log10(args.dm41_max):args.dm41_N*1j,
		np.log10(args.Ue4sq_min):np.log10(args.Ue4sq_max):args.Ue4sq_N*1j,
		np.log10(args.Um4sq_min):np.log10(args.Um4sq_max):args.Um4sq_N*1j,
		np.log10(args.Ut4sq_min):np.log10(args.Ut4sq_max):args.Ut4sq_N*1j,
		]).reshape(4, args.dm41_N*args.Ue4sq_N*args.Um4sq_N*args.Ut4sq_N).T


def write_grid_cfg(args):
	with open("grids/%s/params.cfg"%args.gridname, "w") as _f:
		for a in ["dm41", "Ue4sq", "Um4sq", "Ut4sq"]:
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
	for a in ["dm41", "Ue4sq", "Um4sq", "Ut4sq"]:
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
	return values, fillGrid(Namespace(**values))


def contourplot(xv, yv, values, points,
		fname=None,
		levels=[0., 0.1, 0.3, 0.5, 0.7, 0.9],
		title=None,
		xlab=None, ylab=None,
		xlim=None, ylim=None,
		bfx=-1, bfy=-1, bfup=False,
		lsn_contours=[], lsn_contours_convssq2th=True,
		):
	try:
		zv = points.reshape((len(yv), len(xv)))
	except ValueError:
		print("You probably have to check the params you are using "
			+ "and what you fix, as the shapes of the arrays don't match.\n"
			+ "x:%s, y:%s, z:%s"%(xv.shape, yv.shape, points.shape))
		return
	fig = plt.Figure(figsize=(6,4))
	cf = plt.contourf(xv, yv, zv, levels=levels, cmap=matplotlib.cm.get_cmap('CMRmap'), extend="max")
	cf.cmap.set_over(cmap(0.95))
	for fn in lsn_contours:
		data = np.loadtxt(fn)
		if lsn_contours_convssq2th:
			plt.plot((1.-np.sqrt(1.-data[:,0]))*0.5, data[:,1], '-k')
		else:
			plt.plot(data[:,0], data[:,1], '-k')
	cbar = plt.colorbar(cf)
	if bfx>0 and bfy>0:
		plt.plot(bfx, bfy, color="g", marker=r"$\leftarrow$" if bfup else "*", markersize=10)
	cbar.ax.set_ylabel(r"$\Delta N_{\rm eff}$")
	if title is not None:
		plt.title(title, y=1.02)
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
	plt.tight_layout(rect=(-0.03, -0.05, 1.05, 1.02))
	if fname is not None:
		plt.savefig(fname)
	plt.close()


def call_plot(args):
	def convert_grid(convert_x, convert_y, pts):
		outpts = np.column_stack((pts, np.zeros(np.shape(pts)[0])))
		outpts = np.column_stack((outpts, np.zeros(np.shape(pts)[0])))
		cva = [convert_x, convert_y]
		if "sinsqth14" in cva:
			ix = cva.index("sinsqth14")
			outpts[:, 4+ix] = pts[:,1]
		if "sinsqth24" in cva:
			ix = cva.index("sinsqth24")
			outpts[:, 4+ix] = pts[:,2]/(1-pts[:,1])
		if "sinsqth34" in cva:
			ix = cva.index("sinsqth34")
			outpts[:, 4+ix] = pts[:,3]/(1-pts[:,1]-pts[:,2])
		if "sinsq2th_ee" in cva:
			ix = cva.index("sinsq2th_ee")
			outpts[:, 4+ix] = 4.*pts[:,1]*(1.-pts[:,1])
		if "sinsq2th_mumu" in cva:
			ix = cva.index("sinsq2th_mumu")
			outpts[:, 4+ix] = 4.*pts[:,2]*(1.-pts[:,2])
		if "sinsq2th_emu" in cva:
			ix = cva.index("sinsq2th_emu")
			outpts[:, 4+ix] = 4.*pts[:,1]*pts[:,2]
		return outpts
	mixings, fullgrid, fullobjects = call_read(args)
	fullpoints = list(map(lambda x: safegetattr(x, "Neff", 0.)-args.Neff_active, fullobjects))
	fullpoints = np.asarray(fullpoints)
	cgrid = {}
	if args.par_x == args.par_y:
		print("Cannot plot the same variable in the two axis! %s, %s"%(args.par_x, args.par_y))
		return
	fullgrid = convert_grid(args.par_x, args.par_y, fullgrid)
	for a in indexes.keys():
		if getattr(args, "fix_%s"%a) == "False":
			setattr(args, "fix_%s"%a, False)
		if args.par_x == a or args.par_y == a:
			setattr(args, "fix_%s"%a, False)
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
	if "sinsq" in args.par_x:
		xv = np.unique(smallgrid[:,4])
	else:
		xv = mixings["%s_pts"%args.par_x]
	if "sinsq" in args.par_y:
		xv = np.unique(smallgrid[:,5])
	else:
		yv = mixings["%s_pts"%args.par_y]
	if args.lsn_contours:
		if (args.par_y == "dm41"
				and args.par_x in ["Ue4sq", "sinsqth14", "sinsq2th_ee"]):
			lsn_contours = [
				"/home/gariazzo/data/lsn/globalfit/1801.06469/fig4b-contours/cnt-9973-%d.dat"%i
				for i in [1,2,3]]
			lsn_contours_convssq2th = args.par_x in ["Ue4sq", "sinsqth14"]
	else:
		lsn_contours = []
		lsn_contours_convssq2th = False
	contourplot(
		xv, yv, mixings, smallpoints,
		fname="grids/%s/plots/%s_%s.pdf"%(args.gridname, args.par_x, args.par_y)
			if args.filename == "" else "grids/%s/plots/%s"%(args.gridname, args.filename),
		xlab=labels[args.par_x], ylab=labels[args.par_y],
		title=args.title,
		bfx=args.bestfit_x, bfy=args.bestfit_y, bfup=args.bestfit_upper,
		lsn_contours=lsn_contours,
		lsn_contours_convssq2th=lsn_contours_convssq2th,
		)
	print("\nDone!\n\n")
	return


def call_prepare(args):
	if not os.path.exists("grids/%s/ini/"%args.gridname):
		os.makedirs("grids/%s/ini/"%args.gridname)
	if not os.path.exists("grids/%s/OUT/"%args.gridname):
		os.makedirs("grids/%s/OUT/"%args.gridname)
	if not os.path.exists("grids/%s/plots/"%args.gridname):
		os.makedirs("grids/%s/plots/"%args.gridname)
	files = list(glob.iglob("grids/%s/ini/*.ini"%args.gridname))
	list(map(lambda x: os.remove(x), files))
	for a in ["dm41", "Ue4sq", "Um4sq", "Ut4sq"]:
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
	for dm41, Ue4sq, Um4sq, Ut4sq in grid:
		ssq14 = Ue4sq
		ssq24 = Um4sq/(1.-Ue4sq)
		ssq34 = Ut4sq/(1.-Ue4sq-Um4sq)
		prep = [
			"grids/%s/ini/%s_%s_%s_%s.ini"%(args.gridname, dm41, Ue4sq, Um4sq, Ut4sq),
			"grids/%s/OUT/%s_%s_%s_%s/"%(args.gridname, dm41, Ue4sq, Um4sq, Ut4sq),
			"3+1",
			"damping",
			"--dlsoda_rtol=%s"%args.tolerance,
			"--dlsoda_atol=%s"%args.tolerance,
			"--Nx=%s"%args.Nx,
			"--Ny=%s"%args.Ny,
			"--Nylog=%s"%args.Nylog,
			"--y_cen=1",
			"--x_in=%s"%args.x_in,
			"--default_sterile=None",
			"--dm41=%s"%dm41,
			"--th14=%s"%ssq14,
			"--th24=%s"%ssq24,
			"--th34=%s"%ssq34,
			"--ordering=%s"%args.ordering,
			]
		parser = prepareIni.setParser()
		rargs = parser.parse_args(prep)
		values = prepareIni.getIniValues(rargs)
		prepareIni.writeIni(rargs.inifile, values)
	print("\nTotal number of points: %s"%len(grid))


def call_read(args):
	values, grid = read_grid_cfg(args.gridname)
	objects = []
	missing = 0
	print("\nTotal number of points: %s"%len(grid))
	for dm41, Ue4sq, Um4sq, Ut4sq in grid:
		lab = (r"dm41=%s "%dm41
			+ r"Ue4sq=%s "%Ue4sq
			+ r"Um4sq=%s "%Um4sq
			+ r"Ut4sq=%s "%Ut4sq
			)
		folder = "grids/%s/OUT/%s_%s_%s_%s/"%(args.gridname, dm41, Ue4sq, Um4sq, Ut4sq)
		obj = None
		try:
			obj = NuDensRun(folder, label=lab, nnu=4, rho=False)
		except (IOError, IndexError):
			print("no %s"%lab)
			missing += 1
		else:
			try:
				obj.Neff
			except AttributeError:
				missing += 1
		objects.append(obj)
	print("\nMissing or incomplete points: %s\n"%missing)
	return values, grid, objects


def call_run(args):
	print("submitting the grid %s"%args.gridname)
	if args.failed_only:
		files = []
		for f in list(glob.iglob("grids/%s/OUT/*/z.dat"%args.gridname)):
			if "%s"%args.expected_lines not in subprocess.check_output(['wc', '-l', f]).decode("utf-8"):
				if args.remove_existing:
					shutil.rmtree(f.replace("z.dat", ""), ignore_errors=True)
				files.append(f.replace("OUT", "ini").replace("/z.dat", ".ini"))
	else:
		files = list(glob.iglob("grids/%s/ini/*.ini"%args.gridname))
	current = 0
	for i, f in enumerate(sorted(files)):
		if i >= args.first_index and i < args.last_index:
			current += 1
			os.system(
				"clusterlauncher -N %s_%s -n 1 --openmp -q short-seq -w 6:00:00 bin/nuDens.exe %s"%(
					args.gridname,
					f.split(os.sep)[-1].replace(".ini", ""),
					f,
					))
	print("\nTotal number of runs: %s, submitted: %s"%(len(files), current))


if __name__=='__main__':
	parser = setParser()
	args = parser.parse_args(sys.argv[1:])
	print(args)
	args.func(args)
