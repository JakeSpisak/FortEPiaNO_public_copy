import sys
import os
import argparse
import ast
import glob
import numpy as np
import prepareIni


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
			_f.write("%s: min=%s, max=%s, N=%s\n"%(
				a,
				getattr(args, "%s_min"%a),
				getattr(args, "%s_max"%a),
				getattr(args, "%s_N"%a),
				))


def read_grid_cfg(gridname):
	with open("grids/%s/params.cfg"%args.gridname) as _f:
		lines = _f.read()
		# for a in ["dm41", "ssq14", "ssq24", "ssq34"]:
			# _f.write("%s: min=%s, max=%s, N=%s\n"%(
				# a,
				# getattr(args, "%s_min"%a),
				# getattr(args, "%s_max"%a),
				# getattr(args, "%s_N"%a),
				# ))


def call_plot(args):
	params = readgrid.cfg(args.gridname)
	print("plot")


def call_prepare(args):
	os.makedirs("grids/%s/ini/"%args.gridname, exist_ok=True)
	os.makedirs("grids/%s/OUT/"%args.gridname, exist_ok=True)
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
	params = readgrid.cfg(args.gridname)
	print("read")


def call_run(args):
	files = list(glob.iglob("grids/%s/ini/*.ini"%args.gridname))
	print("submitting the grid %s"%args.gridname)
	for f in files:
		os.system(
			"clusterlauncher -N %s -n 1 --openmp -q short-seq bin/nuDens.exe %s"%(
				f.split(os.sep)[-1].replace(".ini", ""),
				f,
				))


if __name__=='__main__':
	parser = setParser()
	args = parser.parse_args(sys.argv[1:])
	print(args)
	args.func(args)
