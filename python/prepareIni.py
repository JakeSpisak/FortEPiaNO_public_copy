import sys
import os
import argparse

params_active = ["dm21", "th12", "dm31", "th13", "th23"]
params_sterile = ["dm41", "th14", "th24", "th34"]
default_osc_act = {
    "VLC": {
        "no": {
            "dm21": 7.55e-05,
            "th12": 0.32,
            "dm31": 0.00250,
            "th13": 0.0216,
            "th23": 0.547,
        },
        "io": {
            "dm21": 7.55e-05,
            "th12": 0.32,
            "dm31": -0.00242,
            "th13": 0.0222,
            "th23": 0.551,
        },
    },
    "Bari": {
        "no": {
            "dm21": 7.34e-05,
            "dm31": 2.455e-3 + 0.5 * 7.34e-05,
            "th12": 0.304,
            "th13": 0.0214,
            "th23": 0.551,
        },
        "io": {
            "dm21": 7.34e-05,
            "dm31": -2.441e-3 + 0.5 * 7.34e-05,
            "th12": 0.303,
            "th13": 0.0218,
            "th23": 0.557,
        },
    },
    "NuFit": {
        "no": {
            "dm21": 7.39e-05,
            "th12": 0.31,
            "dm31": 0.002525,
            "th13": 0.0224,
            "th23": 0.582,
        },
        "io": {
            "dm21": 7.39e-05,
            "th12": 0.31,
            "dm31": -0.002512 + 7.39e-05,
            "th13": 0.02263,
            "th23": 0.582,
        },
    },
}
default_osc_ster = {
    "Gariazzo&al": {"dm41": 1.29, "th14": 0.01, "th24": 0.0, "th34": 0.0,},
}


def setParser():
    """Prepare the parser for reading the command line arguments

    Output:
        the parser
    """
    parser = argparse.ArgumentParser(prog="prepareIni.py")
    parser.add_argument(
        "inifile",
        metavar="inifilename",
        help="the filename of the ini file where to write the configuration",
    )
    parser.add_argument(
        "outputfolder",
        help="the name of the output folder that will contain the results",
    )
    parser.add_argument(
        "numodel",
        choices=["3nu", "2nu", "3+1", "2+1", "1+1"],
        help="define the neutrino model that must be used",
    )
    parser.add_argument(
        "collisional",
        choices=["zero", "complete", "damping", "diagonal", "yyyw_off", "yyyw_all"],
        help="define the scheme for the collision integrals",
    )
    for c in ["nue", "nunu"]:
        parser.add_argument(
            "--damping_no_%s" % c,
            action="store_true",
            help="disable %s contributions to damping terms" % c,
        )
    parser.add_argument(
        "--qed_corrections",
        choices=["no", "o2", "o3", "o2ln", "o3ln"],
        default="o3",
        help="define which terms must be included "
        + "for the finite-temperature QED corrections "
        + "[O(e^2), O(e^2)+O(e^3) - default, O(e^2)+ln, O(e^2)+O(e^3)ln, or none]",
    )
    parser.add_argument(
        "--ordering",
        choices=["NO", "IO"],
        default="NO",
        help="define the mass ordering for the three active neutrinos "
        + "(not used if you explicitely give the mixing parameters,"
        + " only if you select default values by Valencia, Bari or NuFit global fit)",
    )
    parser.add_argument(
        "--default_active",
        choices=["Bari", "NuFit", "VLC", "None"],
        default="VLC",
        help="define the mixing parameters for the active neutrinos as obtained "
        + "from the Valencia global fit (doi:10.1016/j.physletb.2018.06.019, default), "
        + "Bari group (doi:10.1016/j.ppnp.2018.05.005) or "
        + "NuFit analysis (doi:10.1007/JHEP01(2019)106)",
    )
    parser.add_argument(
        "--default_sterile",
        choices=["Gariazzo&al", "None"],
        default="Gariazzo&al",
        help="define the active-sterile mixing parameters as obtained "
        + "from the Gariazzo et al. global fit (with th24=th34=0)",
    )
    parser.add_argument(
        "--sinsq",
        dest="use_sinsq",
        action="store_false",
        help=r"use the $\sin^2$ of the mixing angles as input",
    )
    parser.add_argument(
        "--dm21", type=float, default=0.0, help=r"define $\Delta m^2_{21}$"
    )
    parser.add_argument(
        "--dm31",
        type=float,
        default=0.0,
        help=r"define $\Delta m^2_{31}$ (pass negative value for inverted ordering)",
    )
    parser.add_argument(
        "--dm41", type=float, default=0.0, help=r"define $\Delta m^2_{41}$"
    )
    parser.add_argument(
        "--th12",
        type=float,
        default=0.0,
        help=r"define $\theta_{12}$ or $\sin^2 \theta_{12}$",
    )
    parser.add_argument(
        "--th13",
        type=float,
        default=0.0,
        help=r"define $\theta_{13}$ or $\sin^2 \theta_{13}$",
    )
    parser.add_argument(
        "--th14",
        type=float,
        default=0.0,
        help=r"define $\theta_{14}$ or $\sin^2 \theta_{14}$",
    )
    parser.add_argument(
        "--th23",
        type=float,
        default=0.0,
        help=r"define $\theta_{23}$ or $\sin^2 \theta_{23}$",
    )
    parser.add_argument(
        "--th24",
        type=float,
        default=0.0,
        help=r"define $\theta_{24}$ or $\sin^2 \theta_{24}$",
    )
    parser.add_argument(
        "--th34",
        type=float,
        default=0.0,
        help=r"define $\theta_{34}$ or $\sin^2 \theta_{34}$",
    )
    parser.add_argument(
        "-V", "--verbose", type=int, default=1, help="define the verbosity of the code"
    )
    parser.add_argument(
        "--verbose_deriv_freq",
        type=int,
        default=100,
        help="print a string stating the current position only after N derivatives",
    )
    parser.add_argument(
        "--Nx", type=int, default=200, help="number of points to save in x"
    )
    parser.add_argument("--x_in", type=float, default=0.001, help="initial value of x")
    parser.add_argument("--x_fin", type=float, default=35, help="final value of x")
    parser.add_argument(
        "--Ny", type=int, default=30, help="number of total points in y"
    )
    parser.add_argument(
        "--Nylog",
        type=int,
        default=5,
        help="number of log-spaced points between y_in and y_cen",
    )
    parser.add_argument("--y_min", type=float, default=0.01, help="minimum value of y")
    parser.add_argument(
        "--y_cen",
        type=float,
        default=1,
        help="value of y where to switch between log- and linear spacing",
    )
    parser.add_argument("--y_max", type=float, default=20, help="maximum value of y")
    parser.add_argument(
        "--dlsoda_atol",
        type=float,
        default=1e-6,
        help="absolute tolerance for all the differential equations in DLSODA. "
        + "See also dlsoda_atol_z, dlsoda_atol_d, dlsoda_atol_o",
    )
    parser.add_argument(
        "--dlsoda_atol_z",
        type=float,
        default=1e-6,
        help="absolute tolerance for the dz/dx, dw/dx differential equations in DLSODA. "
        + "See also dlsoda_atol, dlsoda_atol_d, dlsoda_atol_o",
    )
    parser.add_argument(
        "--dlsoda_atol_d",
        type=float,
        default=1e-6,
        help="absolute tolerance for the differential equations drho_{ii}/dx"
        + "of the diagonal matrix elements in DLSODA. "
        + "See also dlsoda_atol, dlsoda_atol_z, dlsoda_atol_o",
    )
    parser.add_argument(
        "--dlsoda_atol_o",
        type=float,
        default=1e-6,
        help="absolute tolerance for the differential equations drho_{ij}/dx"
        + "of the off-diagonal matrix elements in DLSODA. "
        + "See also dlsoda_atol, dlsoda_atol_z, dlsoda_atol_d",
    )
    parser.add_argument(
        "--dlsoda_rtol", type=float, default=1e-6, help="relative tolerance for DLSODA"
    )
    parser.add_argument(
        "--save_energy_entropy",
        action="store_true",
        help="enable saving the evolution of the energy density and entropy for each component",
    )
    parser.add_argument(
        "--save_fd",
        action="store_true",
        help="enable saving the y grid and the corresponding Fermi-Dirac to fd.dat",
    )
    parser.add_argument(
        "--save_Neff", action="store_true", help="enable saving the evolution of Neff"
    )
    parser.add_argument(
        "--save_nuDens",
        action="store_true",
        help="enable saving the evolution of the full neutrino density matrix",
    )
    parser.add_argument(
        "--save_z",
        action="store_true",
        help="enable saving the evolution of the photon temperature z",
    )
    parser.add_argument(
        "--no_GL",
        action="store_true",
        help="do not use the Gauss-Laguerre method for integrals and for spacing the y points",
    )
    return parser


def oscParams(args):
    """Read part of the settings and prepare the configuration
    for the number of neutrinos and the oscillation parameters

    Parameters:
        args: the output of parse_args

    Output:
        a dictionary
    """
    osc = {}
    osc["use_sinsq"] = "T" if args.use_sinsq else "F"
    if args.numodel in ["3p1", "3+1"]:
        osc["nnu"] = 4
        osc["sterile"] = [False, False, False, True]
        osc["factors"] = [1, 1, 1, 1]
        if args.default_sterile in default_osc_ster.keys():
            for p in params_sterile:
                osc[p] = default_osc_ster[args.default_sterile][p]
        else:
            for p in params_sterile:
                osc[p] = getattr(args, p)
    else:
        for p in params_sterile:
            osc[p] = 0.0
    if args.numodel in ["3p1", "3+1", "3p0", "3+0", "3nu", "3", "2p1", "2+1"]:
        if args.numodel in ["2p1", "2+1"]:
            osc["nnu"] = 3
            osc["sterile"] = [False, False, True]
            osc["factors"] = [1, 2, 1]
        elif args.numodel in ["3p0", "3+0", "3nu", "3"]:
            osc["nnu"] = 3
            osc["sterile"] = [False, False, False]
            osc["factors"] = [1, 1, 1]
        if (
            args.numodel in ["3p1", "3+1", "3p0", "3+0", "3nu", "3",]
            and args.default_active in default_osc_act.keys()
        ):
            for p in params_active:
                osc[p] = default_osc_act[args.default_active][args.ordering.lower()][p]
        else:
            for p in params_active:
                osc[p] = getattr(args, p)
    else:
        osc["dm31"] = 0.0
        osc["th13"] = 0.0
        osc["th23"] = 0.0
    if args.numodel in ["1p1", "1+1", "a+s", "as", "2+0", "2nu", "2"]:
        osc["nnu"] = 2
        if args.numodel in ["a+s", "as", "1p1", "1+1"]:
            osc["sterile"] = [False, True]
            osc["factors"] = [1, 1]
        elif args.numodel in ["2+0", "2nu", "2"]:
            osc["sterile"] = [False, False]
            osc["factors"] = [1, 2]
        osc["dm21"] = args.dm21
        osc["th12"] = args.th12
    return osc


def getIniValues(args):
    """Read the input namespace (or args from argparser)
    and prepare a dictionary for later use
    in the writing of the ini file

    Parameters:
        args: the output of parse_args

    Output:
        a dictionary
    """
    values = oscParams(args)
    for p in [
        "verbose",
        "Nx",
        "x_in",
        "x_fin",
        "Ny",
        "Nylog",
        "y_min",
        "y_cen",
        "y_max",
    ]:
        values[p] = getattr(args, p)
    values["factors"] = "\n".join(
        ["nuFactor%d = %f" % (i + 1, f) for i, f in enumerate(values["factors"])]
    )
    values["sterile"] = "\n".join(
        [
            "sterile%d = %s" % (i + 1, "T" if f else "F")
            for i, f in enumerate(values["sterile"])
        ]
    )
    values["coll_offdiag"] = (
        0
        if args.collisional == "zero"
        else 1
        if args.collisional == "complete"
        else 2
        if args.collisional == "damping"
        else 4
        if args.collisional == "yyyw_off"
        else 5
        if args.collisional == "yyyw_all"
        else 3
    )
    if args.qed_corrections == "no":
        values["ftqed_temperature_corr"] = "F"
        values["ftqed_ord3"] = "F"
        values["ftqed_log_term"] = "F"
    else:
        values["ftqed_temperature_corr"] = "T"
        values["ftqed_ord3"] = "T" if "o3" in args.qed_corrections else "F"
        values["ftqed_log_term"] = "T" if "ln" in args.qed_corrections else "F"
    if any(
        [
            a != 1e-6
            for a in [args.dlsoda_atol_z, args.dlsoda_atol_d, args.dlsoda_atol_o]
        ]
    ):
        values["dlsoda_atol"] = (
            "dlsoda_atol_z = %s\n" % args.dlsoda_atol_z
            + "dlsoda_atol_d = %s\n" % args.dlsoda_atol_d
            + "dlsoda_atol_o = %s\n" % args.dlsoda_atol_o
        )
    else:
        values["dlsoda_atol"] = (
            "dlsoda_atol_z = %s\n" % args.dlsoda_atol
            + "dlsoda_atol_d = %s\n" % args.dlsoda_atol
            + "dlsoda_atol_o = %s\n" % args.dlsoda_atol
        )
    values["dlsoda_rtol"] = args.dlsoda_rtol
    values["folder"] = args.outputfolder
    values["Nprintderivs"] = args.verbose_deriv_freq
    for p in [
        "damping_no_nue",
        "damping_no_nunu",
        "save_energy_entropy",
        "save_fd",
        "save_Neff",
        "save_nuDens",
        "save_z",
    ]:
        values[p] = "T" if getattr(args, p) else "F"
    values["use_GL"] = "F" if args.no_GL else "T"
    return values


def writeIni(filename, values):
    """Use the information already prepared
    to write an ini file for FortEPiaNO

    Parameters:
        filename: the name of the output ini file
        values: a dictionary with the settings for the ini file
    """
    iniText = """###run setttings
flavorNumber = {nnu:}

{factors:}
{sterile:}

givesinsq = {use_sinsq:}
theta12= {th12:}
dm21 = {dm21:}
theta13 = {th13:}
theta23 = {th23:}
dm31 = {dm31:}
theta14 = {th14:}
theta24 = {th24:}
theta34 = {th34:}
dm41 = {dm41:}

collision_offdiag = {coll_offdiag:}
ftqed_temperature_corr = {ftqed_temperature_corr:}
ftqed_ord3 = {ftqed_ord3:}
ftqed_log_term = {ftqed_log_term:}

Nx = {Nx:}
x_in = {x_in:}
x_fin = {x_fin:}

use_gauss_laguerre = {use_GL:}
Ny = {Ny:}
Nylog = {Nylog:}
y_min = {y_min:}
y_cen = {y_cen:}
y_max = {y_max:}

outputFolder = {folder:}
checkpoint = T

save_fd = {save_fd:}
save_Neff = {save_Neff:}
save_nuDens_evolution = {save_nuDens:}
save_z_evolution = {save_z:}
save_energy_entropy_evolution = {save_energy_entropy:}

{dlsoda_atol:}
dlsoda_rtol = {dlsoda_rtol:}

verbose = {verbose:}
Nprintderivs = {Nprintderivs:}
""".format(
        **values
    )
    print("Writing to %s" % filename)
    dirname = os.path.dirname(os.path.abspath(filename))
    if not os.path.exists(dirname):
        os.mkdir(dirname)
    with open(filename, "w") as _f:
        _f.write(iniText)


if __name__ == "__main__":
    parser = setParser()
    args = parser.parse_args(sys.argv[1:])
    values = getIniValues(args)
    writeIni(args.inifile, values)
