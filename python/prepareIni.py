import sys
import argparse


def setParser():
	parser = argparse.ArgumentParser(prog='prepareIni.py')
	parser.add_argument(
		'inifile',
		metavar='inifilename',
		help='the filename of the ini file where to write the configuration'
		)
	parser.add_argument(
		'outputfolder',
		help='the name of the output folder that will contain the results'
		)
	parser.add_argument(
		'numodel',
		choices=["3nu", "2nu", "3+1", "2+1", "1+1"],
		help='define the neutrino model that must be used'
		)
	parser.add_argument(
		'collisional',
		choices=["zero", "complete", "damping", "diagonal"],
		help='define the scheme for the collision integrals'
		)
	parser.add_argument(
		'--ordering',
		choices=["NO", "IO"],
		default="NO",
		help='define the mass ordering for the three active neutrinos'
		)
	parser.add_argument(
		'--default_active',
		nargs=1,
		choices=["VLC", ""],
		default="VLC",
		help='define the mixing parameters for the active neutrinos as obtained from the Valencia global fit'
		)
	parser.add_argument(
		'--default_sterile',
		nargs=1,
		choices=["Gariazzo&al", ""],
		default="Gariazzo&al",
		help='define the active-sterile mixing parameters as obtained from the Gariazzo et al. global fit (with th24=th34=0)'
		)
	parser.add_argument(
		'--sinsq',
		dest='use_sinsq',
		action='store_false',
		help=r'use the $\sin^2$ of the mixing angles as input'
		)
	parser.add_argument(
		'--dm21',
		type=float,
		default=0.,
		help=r'define $\Delta m^2_{21}$'
		)
	parser.add_argument(
		'--dm31',
		type=float,
		default=0.,
		help=r'define $|\Delta m^2_{31}|$'
		)
	parser.add_argument(
		'--dm41',
		type=float,
		default=0.,
		help=r'define $\Delta m^2_{41}$'
		)
	parser.add_argument(
		'--th12',
		type=float,
		default=0.,
		help=r'define $\theta_{12}$ or $\sin^2 \theta_{12}$'
		)
	parser.add_argument(
		'--th13',
		type=float,
		default=0.,
		help=r'define $\theta_{13}$ or $\sin^2 \theta_{13}$'
		)
	parser.add_argument(
		'--th14',
		type=float,
		default=0.,
		help=r'define $\theta_{14}$ or $\sin^2 \theta_{14}$'
		)
	parser.add_argument(
		'--th23',
		type=float,
		default=0.,
		help=r'define $\theta_{23}$ or $\sin^2 \theta_{23}$'
		)
	parser.add_argument(
		'--th24',
		type=float,
		default=0.,
		help=r'define $\theta_{24}$ or $\sin^2 \theta_{24}$'
		)
	parser.add_argument(
		'--th34',
		type=float,
		default=0.,
		help=r'define $\theta_{34}$ or $\sin^2 \theta_{34}$'
		)
	parser.add_argument(
		'-V',
		'--verbose',
		type=int,
		default=1,
		help='define the verbosity of the code',
		)
	parser.add_argument(
		'--verbose_deriv_freq',
		type=int,
		default=100,
		help='print a string stating the current position only after N derivatives',
		)
	parser.add_argument(
		'--Nx',
		type=int,
		default=200,
		help='number of points to save in x',
		)
	parser.add_argument(
		'--x_in',
		type=float,
		default=0.001,
		help='initial value of x',
		)
	parser.add_argument(
		'--x_fin',
		type=float,
		default=35,
		help='final value of x',
		)
	parser.add_argument(
		'--Ny',
		type=int,
		default=40,
		help='number of total points in y',
		)
	parser.add_argument(
		'--Nylog',
		type=int,
		default=5,
		help='number of log-spaced points between y_in and y_cen',
		)
	parser.add_argument(
		'--y_min',
		type=float,
		default=0.01,
		help='minimum value of y',
		)
	parser.add_argument(
		'--y_cen',
		type=float,
		default=1,
		help='value of y where to switch between log- and linear spacing',
		)
	parser.add_argument(
		'--y_max',
		type=float,
		default=20,
		help='maximum value of y',
		)
	parser.add_argument(
		'--dlsoda_atol',
		type=float,
		default=1e-6,
		help='absolute tolerance for DLSODA',
		)
	parser.add_argument(
		'--dlsoda_rtol',
		type=float,
		default=1e-6,
		help='relative tolerance for DLSODA',
		)
	return parser


def oscParams(args):
	osc = {}
	osc["use_sinsq"] = "T" if args.use_sinsq else "F"
	osc["ordering"] = "T"
	if args.numodel in ["3p1", "3+1"]:
		osc["nnu"] = 4
		osc["sterile"] = [False, False, False, True]
		osc["factors"] = [1, 1, 1, 1]
		if args.default_sterile == "Gariazzo&al":
			osc["dm41"] = 1.29
			osc["th14"] = 0.01
			osc["th24"] = 0.
			osc["th34"] = 0.
		else:
			osc["dm41"] = args.dm41
			osc["th14"] = args.th14
			osc["th24"] = args.th24
			osc["th34"] = args.th34
	else:
		osc["dm41"] = 0.
		osc["th14"] = 0.
		osc["th24"] = 0.
		osc["th34"] = 0.
	if args.numodel in ["3p1", "3+1", "3p0", "3+0", "3nu", "3", "2p1", "2+1"]:
		if args.numodel in ["2p1", "2+1"]:
			osc["nnu"] = 3
			osc["sterile"] = [False, False, True]
			osc["factors"] = [1, 2, 1]
		elif args.numodel in ["3p0", "3+0", "3nu", "3"]:
			osc["nnu"] = 3
			osc["sterile"] = [False, False, False]
			osc["factors"] = [1, 1, 1]
		if args.numodel in ["3p1", "3+1", "3p0", "3+0", "3nu", "3"] \
				and args.default_active == "VLC":
			if args.ordering.lower() in ["no", "nh", "normal"]:
				osc["ordering"] = "T"
				osc["dm21"] = 0.32
				osc["dm31"] = 0.00250
				osc["th12"] = 7.55e-05
				osc["th13"] = 0.0216
				osc["th23"] = 0.547
			else:
				osc["ordering"] = "F"
				osc["dm21"] = 0.32
				osc["dm31"] = 0.00242
				osc["th12"] = 7.55e-05
				osc["th13"] = 0.0222
				osc["th23"] = 0.551
		else:
			if args.ordering.lower() in ["no", "nh", "normal"]:
				osc["ordering"] = "T"
			else:
				osc["ordering"] = "F"
			osc["dm21"] = args.dm21
			osc["dm31"] = args.dm31
			osc["th12"] = args.th12
			osc["th13"] = args.th13
			osc["th23"] = args.th23
	else:
		osc["dm31"] = 0.
		osc["th13"] = 0.
		osc["th23"] = 0.
	if args.numodel in ["1p1", "1+1", "a+s", "as", "2+0", "2nu", "2"]:
		osc["nnu"] = 2
		if args.numodel in ["1p1", "1+1"]:
			osc["sterile"] = [False, True]
			osc["factors"] = [3, 1]
		elif args.numodel in ["a+s", "as", "2+0", "2nu", "2"]:
			osc["sterile"] = [False, False]
			osc["factors"] = [1, 2]
		osc["dm21"] = args.dm21
		osc["th12"] = args.th12
	return osc


def getIniValues(args):
	values = oscParams(args)
	values["verbose"] = args.verbose
	values["factors"] = "\n".join(
		["nuFactor%d = %f"%(i+1, f) for i, f in enumerate(values["factors"])])
	values["sterile"] = "\n".join(
		["sterile%d = %s"%(i+1, "T" if f else "F")
			for i, f in enumerate(values["sterile"])])
	values["coll_offdiag"] = \
		0 if args.collisional == "zero" \
		else 1 if args.collisional == "complete" \
		else 2 if args.collisional == "damping" \
		else 3
	values["Nx"] = args.Nx
	values["x_in"] = args.x_in
	values["x_fin"] = args.x_fin
	values["Ny"] = args.Ny
	values["Nylog"] = args.Nylog
	values["y_min"] = args.y_min
	values["y_cen"] = args.y_cen
	values["y_max"] = args.y_max
	values["dlsoda_atol"] = args.dlsoda_atol
	values["dlsoda_rtol"] = args.dlsoda_rtol
	values["folder"] = args.outputfolder
	values["Nprintderivs"] = args.verbose_deriv_freq
	return values


def writeIni(filename, values):
	iniText = """###run setttings
flavorNumber = {nnu:}

{factors:}
{sterile:}

massOrdering = {ordering:}

givesinsq = {use_sinsq:}
theta12= {th12:}
dm12 = {dm21:}
theta13 = {th13:}
theta23 = {th23:}
dm31 = {dm31:}
theta14 = {th14:}
theta24 = {th24:}
theta34 = {th34:}
dm14 = {dm41:}

collision_offdiag = {coll_offdiag:}
dme2_temperature_corr = T

Nx = {Nx:}
x_in = {x_in:}
x_fin = {x_fin:}

Ny = {Ny:}
Nylog = {Nylog:}
y_min = {y_min:}
y_cen = {y_cen:}
y_max = {y_max:}

checkpoint = T
verbose = {verbose:}
Nprintderivs = {Nprintderivs:}
dlsoda_atol = {dlsoda_atol:}
dlsoda_rtol = {dlsoda_rtol:}
outputFolder = {folder:}
""".format(**values)
	print("Writing to %s"%filename)
	with open(filename, "w") as _f:
		_f.write(iniText)


if __name__=='__main__':
	parser = setParser()
	args = parser.parse_args(sys.argv[1:])
	values = getIniValues(args)
	writeIni(args.inifile, values)