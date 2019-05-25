import sys
import os
import argparse
import ast
import glob
import re
import numpy as np
import shutil
import subprocess

try:
    from scipy.interpolate import interpn
    import matplotlib

    matplotlib.use("agg")
    import matplotlib.pyplot as plt
except ImportError:
    print("missing modules, many things will not work")
else:
    cmap = matplotlib.cm.get_cmap("CMRmap")
try:
    import ternary
    import ternary.helpers
except ImportError:
    print("ternary is not available, some functions are disabled")
    tern = False
    ternary = None
else:
    tern = True
import prepareIni

try:
    from fortepianoOutput import (
        colors,
        markers,
        styles,
        finalizePlot,
        stripRepeated,
        FortEPiaNORun,
    )
except ImportError:
    print("missing modules, many things will not work")


Neff_default = 3.043
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


def sortFiles(files, from_heaviest):
    return (
        sorted(
            files,
            key=lambda x: [
                float(f) for f in x.replace(".ini", "").split("/")[-1].split("_")
            ],
            reverse=True,
        )
        if from_heaviest
        else sorted(files)
    )


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
    parser = argparse.ArgumentParser(prog="dogrid3p1.py")
    parser.add_argument("gridname", help="the name of the grid you want to use")
    parser.add_argument(
        "--verbose",
        action="store_true",
        help="print information on the grid parameters, on the points that are missing and so on",
    )
    subparsers = parser.add_subparsers(help="sub-command help", dest="cmd")
    subparsers.required = True

    parser_fill = subparsers.add_parser(
        "fill",
        help="fill points in the grid that do not need to be computed (you already know that Neff=4.05)",
    )
    parser_fill.add_argument(
        "--fill_up_val",
        type=float,
        default=4.05,
        help="reference value of Neff when fully thermalized 3+1 is considered",
    )
    parser_fill.add_argument(
        "--fill_down_val",
        type=float,
        default=3.044,
        help="reference value of Neff when only active neutrinos are considered",
    )
    parser_fill.add_argument(
        "--fill_up_thres",
        type=float,
        default=4.0,
        help="threshold value of Neff for filling above when fully thermalized 3+1 is considered",
    )
    parser_fill.add_argument(
        "--fill_down_thres",
        type=float,
        default=3.1,
        help="threshold value of Neff for filling below when only active neutrinos are considered",
    )
    parser_fill.set_defaults(func=call_fill)

    parser_plot = subparsers.add_parser("plot", help="plot some output")
    parser_plot.add_argument(
        "--Neff_ref",
        type=float,
        default=0.0,
        help="reference value of Neff when only active neutrinos are considered",
    )
    parser_plot.add_argument(
        "--par_x",
        choices=[
            "Ue4sq",
            "Um4sq",
            "Ut4sq",
            "sinsqth14",
            "sinsqth24",
            "sinsqth34",
            "sinsq2th_ee",
            "sinsq2th_mumu",
            "sinsq2th_emu",
        ],
        default="Um4sq",
        help="parameter for the x axis",
    )
    parser_plot.add_argument(
        "--par_y",
        choices=[
            "dm41",
            "Ue4sq",
            "Um4sq",
            "Ut4sq",
            "sinsqth14",
            "sinsqth24",
            "sinsqth34",
            "sinsq2th_ee",
            "sinsq2th_mumu",
            "sinsq2th_emu",
        ],
        default="dm41",
        help="parameter for the y axis",
    )
    for a, v in [["dm41", False], ["Ue4sq", 0], ["Um4sq", False], ["Ut4sq", 0]]:
        parser_plot.add_argument(
            "--fix_%s" % a,
            default=v,
            help="when requesting a subplot, fix %s to the given index (not value!)"
            % a,
        )
    parser_plot.add_argument("--title", default="", help="title for the plot")
    parser_plot.add_argument(
        "--textbox", default="", help="text to include as a textbox in the figure"
    )
    parser_plot.add_argument(
        "--filename", default="", help="name of the file where to save the plot"
    )
    parser_plot.add_argument(
        "--bestfit_x",
        type=float,
        default=-1,
        help="position of the best-fit in x. Set to negative to ignore",
    )
    parser_plot.add_argument(
        "--bestfit_y",
        type=float,
        default=-1,
        help="position of the best-fit in y. Set to negative to ignore",
    )
    parser_plot.add_argument(
        "--bestfit_upper",
        action="store_true",
        help="The angle specified in the bestfit is only an upper limit",
    )
    parser_plot.add_argument(
        "--lsn_contours", action="store_true", help="Plot contours from global fit"
    )
    parser_plot.add_argument(
        "--colorbar", action="store_false", help="Disable plot contours from global fit"
    )
    parser_plot.add_argument(
        "--colorbar_fname",
        default="",
        help="Name for the file where to save a separate colorbar",
    )
    parser_plot.add_argument(
        "--xmin",
        type=float,
        default=-1,
        help="lower limit for the x axis. Set to negative to ignore",
    )
    parser_plot.add_argument(
        "--xmax",
        type=float,
        default=-1,
        help="upper limit for the x axis. Set to negative to ignore",
    )
    parser_plot.add_argument(
        "--ymin",
        type=float,
        default=-1,
        help="lower limit for the y axis. Set to negative to ignore",
    )
    parser_plot.add_argument(
        "--ymax",
        type=float,
        default=-1,
        help="upper limit for the y axis. Set to negative to ignore",
    )
    parser_plot.set_defaults(func=call_plot)

    parser_prepare = subparsers.add_parser(
        "set", help="set the grid and create ini files"
    )
    parser_prepare.add_argument(
        "--Nx",
        type=int,
        default=1000,
        help="number of points in x where to save the output",
    )
    parser_prepare.add_argument(
        "--no_GL",
        action="store_true",
        help="do not use the Gauss-Laguerre method for integrals and for spacing the y points",
    )
    parser_prepare.add_argument("--Ny", type=int, default=20, help="number of momenta")
    parser_prepare.add_argument(
        "--y_cen",
        type=float,
        default=1,
        help="value of y where to switch between log and linear scale",
    )
    parser_prepare.add_argument(
        "--Nylog",
        type=int,
        default=5,
        help="number of log-spaced points in 0.01 < y < 1",
    )
    parser_prepare.add_argument(
        "--save_energy_entropy",
        action="store_true",
        help="enable saving the evolution of the energy density and entropy for each component",
    )
    parser_prepare.add_argument(
        "--save_fd",
        action="store_true",
        help="enable saving the y grid and the corresponding Fermi-Dirac to fd.dat",
    )
    parser_prepare.add_argument(
        "--save_Neff", action="store_true", help="enable saving the evolution of Neff"
    )
    parser_prepare.add_argument(
        "--save_nuDens",
        action="store_true",
        help="enable saving the evolution of the full neutrino density matrix",
    )
    parser_prepare.add_argument(
        "--save_z",
        action="store_true",
        help="enable saving the evolution of the photon temperature z",
    )
    parser_prepare.add_argument(
        "--x_in", type=float, default=8e-4, help="initial value in x"
    )
    parser_prepare.add_argument(
        "--ternary",
        action="store_true",
        help="activate if you want the grid for ternary plots on the angles [will save a lot of calculations, as it scales with (N+2)(N+1)/2 instead of N^3 for the Ui4sq]",
    )
    parser_prepare.add_argument(
        "--addExp",
        type=float,
        default=0,
        help="Add this number to the exponent, only for ternary grid",
    )
    parser_prepare.add_argument(
        "--tolerance", type=float, default=1e-5, help="tolerance for DLSODA"
    )
    for a, l, mi, ma, n in [
        ["dm41", r"\Delta m^2_{41}", 1e-5, 1e2, 8],
        ["Ue4sq", r"|U_{e4}|^2", 1e-6, 0.1, 6],
        ["Um4sq", r"|U_{\mu4}|^2", 1e-6, 0.1, 6],
        ["Ut4sq", r"|U_{\tau4}|^2", 1e-6, 0.1, 6],
    ]:
        parser_prepare.add_argument(
            "--%s_min" % a, type=float, default=mi, help=r"minimum value of %s" % a
        )
        parser_prepare.add_argument(
            "--%s_max" % a, type=float, default=ma, help=r"maximum value of %s" % a
        )
        parser_prepare.add_argument(
            "--%s_N" % a, type=int, default=n, help=r"number of points in %s" % a
        )
    parser_prepare.add_argument(
        "--ordering",
        choices=["NO", "IO"],
        default="NO",
        help="define the mass ordering for the three active neutrinos",
    )
    parser_prepare.add_argument(
        "--rename", action="store_true", help="rename old to new grid point format"
    )
    parser_prepare.set_defaults(func=call_prepare)

    parser_read = subparsers.add_parser(
        "read", help="read the output files and print resume"
    )
    parser_read.add_argument(
        "--write",
        action="store_true",
        help="save the Neff values in the grid points to a file",
    )
    parser_read.set_defaults(func=call_read)

    parser_run = subparsers.add_parser("run", help="submit the jobs in the grid")
    parser_run.add_argument(
        "-f",
        "--failed_only",
        action="store_true",
        help="resubmit only failed or incomplete runs",
    )
    parser_run.add_argument(
        "-l", "--local", action="store_true", help="run the jobs locally with os.system"
    )
    parser_run.add_argument(
        "-q",
        "--queue",
        default="short-seq",
        help="select the name of the queue where to submit the runs",
    )
    parser_run.add_argument(
        "-r",
        "--remove_existing",
        action="store_true",
        help="remove the folder before resubmitting the failed run",
    )
    parser_run.add_argument(
        "-s",
        "--first-index",
        type=int,
        default=0,
        help="starting index for the runs to submit (if there are a lot of runs)",
    )
    parser_run.add_argument(
        "-t",
        "--last-index",
        type=int,
        default=1000,
        help="ending index for the runs to submit (if there are a lot of runs)",
    )
    parser_run.add_argument(
        "-w",
        "--walltime_hours",
        type=int,
        default=2,
        help="maximum number of hours before killing the job",
    )
    parser_run.add_argument(
        "-m",
        "--walltime_minutes",
        type=int,
        default=0,
        help="maximum number of minutes before killing the job",
    )
    parser_run.add_argument(
        "--from_heaviest",
        action="store_true",
        help="submit the runs starting from the higher masses",
    )
    parser_run.set_defaults(func=call_run)

    parser_ternary = subparsers.add_parser("ternary", help="do a ternary plot")
    parser_ternary.add_argument(
        "--Neff_ref",
        type=float,
        default=0.0,
        help="reference value of Neff when only active neutrinos are considered",
    )
    parser_ternary.add_argument(
        "--fix_dm41",
        type=int,
        default=0,
        help="fix dm41 to the given index (not value!)",
    )
    parser_ternary.add_argument(
        "--filename", default="", help="name of the file where to save the plot"
    )
    parser_ternary.add_argument(
        "--style",
        choices=["heatmap", "scatter"],
        default="heatmap",
        help="which type of plot to use to show Delta Neff",
    )
    parser_ternary.set_defaults(func=call_ternary)
    return parser


def fillGrid(args):
    if tern and hasattr(args, "ternary") and args.ternary:
        pts = []
        try:
            a = args.addExp
        except AttributeError:
            a = 0
        minv = np.log10(args.Ue4sq_min)
        maxv = np.log10(args.Ue4sq_max)
        basescale = maxv - minv
        n = args.Ue4sq_N
        N = int(n / basescale)
        frac = 1.0 / N
        for dm41 in np.logspace(
            np.log10(args.dm41_min), np.log10(args.dm41_max), args.dm41_N
        ):
            for pt in ternary.helpers.simplex_iterator(n):
                pts.append(
                    [
                        dm41,
                        10 ** (a + frac * (pt[0] - n - N)),
                        10 ** (a + frac * (pt[1] - n - N)),
                        10 ** (a + frac * (pt[2] - n - N)),
                    ]
                )
        return np.asarray(pts)
    else:
        return (
            10
            ** (
                np.mgrid[
                    np.log10(args.dm41_min) : np.log10(args.dm41_max) : args.dm41_N
                    * 1j,
                    np.log10(args.Ue4sq_min) : np.log10(args.Ue4sq_max) : args.Ue4sq_N
                    * 1j,
                    np.log10(args.Um4sq_min) : np.log10(args.Um4sq_max) : args.Um4sq_N
                    * 1j,
                    np.log10(args.Ut4sq_min) : np.log10(args.Ut4sq_max) : args.Ut4sq_N
                    * 1j,
                ]
            )
            .reshape(4, args.dm41_N * args.Ue4sq_N * args.Um4sq_N * args.Ut4sq_N)
            .T
        )


def write_grid_cfg(args):
    with open("grids/%s/params.cfg" % args.gridname, "w") as _f:
        for a in ["dm41", "Ue4sq", "Um4sq", "Ut4sq"]:
            _f.write(
                "%s: min=%s, max=%s, N=%s\n"
                % (
                    a,
                    getattr(args, "%s_min" % a),
                    getattr(args, "%s_max" % a),
                    getattr(args, "%s_N" % a),
                )
            )
        if hasattr(args, "ternary") and args.ternary:
            _f.write("ternary\n")


def read_grid_cfg(gridname):
    with open("grids/%s/params.cfg" % gridname) as _f:
        text = _f.readlines()
    values = {"ternary": False}
    for l in text:
        if "ternary" in l:
            values["ternary"] = True
    for a in ["dm41", "Ue4sq", "Um4sq", "Ut4sq"]:
        for l in text:
            res = re.match(
                "%s: min=([0-9\.e\+\-]+), max=([0-9\.e\+\-]+), N=([0-9]+)\n" % a, l
            )
            if res:
                try:
                    values["%s_min" % a] = float(res.group(1))
                    values["%s_max" % a] = float(res.group(2))
                    values["%s_N" % a] = int(res.group(3))
                except (IndexError, TypeError, ValueError):
                    print("cannot read line! %s" % l)
                else:
                    if values["%s_min" % a] == 0 or values["%s_max" % a] == 0:
                        values["%s_pts" % a] = np.linspace(
                            values["%s_min" % a],
                            values["%s_max" % a],
                            values["%s_N" % a]
                            + (1 if (values["ternary"] and a != "dm41") else 0),
                        )
                    else:
                        values["%s_pts" % a] = np.logspace(
                            np.log10(values["%s_min" % a]),
                            np.log10(values["%s_max" % a]),
                            values["%s_N" % a]
                            + (1 if (values["ternary"] and a != "dm41") else 0),
                        )
    return values, fillGrid(Namespace(**values))


def contourplot(
    xv,
    yv,
    values,
    points,
    fname=None,
    levels=[3.0, 3.1, 3.3, 3.5, 3.7, 3.9],
    title=None,
    xlab=None,
    ylab=None,
    xlim=None,
    ylim=None,
    bfx=-1,
    bfy=-1,
    bfup=False,
    lsn_contours=[],
    lsn_contours_convssq2th=True,
    colorbar=True,
    colorbar_fname=None,
    textbox=None,
    cbarlabel=r"$N_{\rm eff}$",
    tightrect=(-0.03, -0.05, 1.05, 1.02),
):
    try:
        zv = points.reshape((len(yv), len(xv)))
    except ValueError:
        print(
            "You probably have to check the params you are using "
            + "and what you fix, as the shapes of the arrays don't match.\n"
            + "x:%s, y:%s, z:%s" % (xv.shape, yv.shape, points.shape)
        )
        return
    fig = plt.figure(figsize=(4.5, 3.0))
    cf = plt.contourf(
        xv, yv, zv, levels=levels, cmap=matplotlib.cm.get_cmap("CMRmap"), extend="both"
    )
    cf.cmap.set_over(cmap(0.95))
    cf.cmap.set_under(cmap(0.0))
    for fn in lsn_contours:
        data = np.loadtxt(fn)
        if lsn_contours_convssq2th:
            plt.plot((1.0 - np.sqrt(1.0 - data[:, 0])) * 0.5, data[:, 1], "-k")
        else:
            plt.plot(data[:, 0], data[:, 1], "-k")
    cbar = plt.colorbar(cf)
    cbar.ax.set_ylabel(cbarlabel)
    if bfx > 0 and bfy > 0:
        plt.plot(
            bfx, bfy, color="g", marker=r"$\leftarrow$" if bfup else "*", markersize=10
        )
    if title is not None:
        plt.title(title, y=1.02)
    if textbox is not None:
        plt.text(
            8e-2,
            40,
            textbox,
            color="k",
            horizontalalignment="right",
            bbox=dict(facecolor="white", alpha=0.6, boxstyle="round"),
        )
    ax = plt.gca()
    ax.tick_params(
        "both",
        which="both",
        direction="out",
        left=True,
        right=True,
        top=True,
        bottom=True,
    )
    ax.tick_params("both", which="major", direction="out", length=8)
    plt.xscale("log")
    plt.yscale("log")
    if xlab is not None:
        plt.xlabel(xlab, fontsize="large")
    if ylab is not None:
        plt.ylabel(ylab, fontsize="large")
    if xlim is not None:
        plt.xlim(xlim)
    if ylim is not None:
        plt.ylim(ylim)
    if not colorbar:
        cbar.remove()
        plt.draw()
        plt.tight_layout(rect=(-0.04, -0.07, 1.19, 1.03))
    else:
        plt.tight_layout(rect=tightrect)
    if fname is not None:
        plt.savefig(fname)
    if colorbar_fname is not None:
        ax.set_visible(False)
        plt.close()
        plt.figure(figsize=(6.0, 1.0))
        cbar = plt.colorbar(cf, orientation="horizontal")
        cbar.ax.set_xlabel(cbarlabel, fontsize="xx-large")
        plt.tight_layout(rect=(-0.02, -0.08, 1.02, 3.0))
        plt.savefig(colorbar_fname)
    plt.close()


def call_fill(args):
    values, grid = read_grid_cfg(args.gridname)
    nsvals = Namespace(**values)
    objects = []
    missing = 0
    filled = 0
    if args.verbose:
        print("\nTotal number of points: %s" % len(grid))
    # first, read existing points
    Neffpoints = {}
    curr = {}
    for currpt in grid:
        curr["dm41"] = currpt[0]
        curr["Ue4sq"] = currpt[1]
        curr["Um4sq"] = currpt[2]
        curr["Ut4sq"] = currpt[3]
        lab = (
            r"dm41=%.5e " % curr["dm41"]
            + r"Ue4sq=%.5e " % curr["Ue4sq"]
            + r"Um4sq=%.5e " % curr["Um4sq"]
            + r"Ut4sq=%.5e " % curr["Ut4sq"]
        )
        string = "%.5e_%.5e_%.5e_%.5e" % (
            curr["dm41"],
            curr["Ue4sq"],
            curr["Um4sq"],
            curr["Ut4sq"],
        )
        folder = "grids/%s/OUT/%s/" % (args.gridname, string)
        try:
            obj = FortEPiaNORun(
                folder, label=lab, nnu=4, rho=False, verbose=args.verbose
            )
        except (IOError, IndexError):
            if args.verbose:
                print("no %s" % lab)
        else:
            try:
                Neffpoints[string] = obj.Neff
            except AttributeError:
                Neffpoints[string] = None
    # define grids over which to scan
    order = [
        ["dm41", "Ue4sq", "Um4sq", "Ut4sq"],
        ["Ut4sq", "dm41", "Ue4sq", "Um4sq"],
        ["Um4sq", "Ut4sq", "dm41", "Ue4sq"],
        ["Ue4sq", "Um4sq", "Ut4sq", "dm41"],
    ]
    grids = [
        10
        ** (
            np.mgrid[
                np.log10(nsvals.dm41_min) : np.log10(nsvals.dm41_max) : nsvals.dm41_N
                * 1j,
                np.log10(nsvals.Ue4sq_min) : np.log10(nsvals.Ue4sq_max) : nsvals.Ue4sq_N
                * 1j,
                np.log10(nsvals.Um4sq_min) : np.log10(nsvals.Um4sq_max) : nsvals.Um4sq_N
                * 1j,
                np.log10(nsvals.Ut4sq_min) : np.log10(nsvals.Ut4sq_max) : nsvals.Ut4sq_N
                * 1j,
            ]
        )
        .reshape(4, nsvals.dm41_N * nsvals.Ue4sq_N * nsvals.Um4sq_N * nsvals.Ut4sq_N)
        .T,
        10
        ** (
            np.mgrid[
                np.log10(nsvals.Ut4sq_min) : np.log10(nsvals.Ut4sq_max) : nsvals.Ut4sq_N
                * 1j,
                np.log10(nsvals.dm41_min) : np.log10(nsvals.dm41_max) : nsvals.dm41_N
                * 1j,
                np.log10(nsvals.Ue4sq_min) : np.log10(nsvals.Ue4sq_max) : nsvals.Ue4sq_N
                * 1j,
                np.log10(nsvals.Um4sq_min) : np.log10(nsvals.Um4sq_max) : nsvals.Um4sq_N
                * 1j,
            ]
        )
        .reshape(4, nsvals.dm41_N * nsvals.Ue4sq_N * nsvals.Um4sq_N * nsvals.Ut4sq_N)
        .T,
        10
        ** (
            np.mgrid[
                np.log10(nsvals.Um4sq_min) : np.log10(nsvals.Um4sq_max) : nsvals.Um4sq_N
                * 1j,
                np.log10(nsvals.Ut4sq_min) : np.log10(nsvals.Ut4sq_max) : nsvals.Ut4sq_N
                * 1j,
                np.log10(nsvals.dm41_min) : np.log10(nsvals.dm41_max) : nsvals.dm41_N
                * 1j,
                np.log10(nsvals.Ue4sq_min) : np.log10(nsvals.Ue4sq_max) : nsvals.Ue4sq_N
                * 1j,
            ]
        )
        .reshape(4, nsvals.dm41_N * nsvals.Ue4sq_N * nsvals.Um4sq_N * nsvals.Ut4sq_N)
        .T,
        10
        ** (
            np.mgrid[
                np.log10(nsvals.Ue4sq_min) : np.log10(nsvals.Ue4sq_max) : nsvals.Ue4sq_N
                * 1j,
                np.log10(nsvals.Um4sq_min) : np.log10(nsvals.Um4sq_max) : nsvals.Um4sq_N
                * 1j,
                np.log10(nsvals.Ut4sq_min) : np.log10(nsvals.Ut4sq_max) : nsvals.Ut4sq_N
                * 1j,
                np.log10(nsvals.dm41_min) : np.log10(nsvals.dm41_max) : nsvals.dm41_N
                * 1j,
            ]
        )
        .reshape(4, nsvals.dm41_N * nsvals.Ue4sq_N * nsvals.Um4sq_N * nsvals.Ut4sq_N)
        .T,
    ]
    # start checking the grid content
    for fill_up in [False, True]:
        for j, curr_ord in enumerate(order):
            print(curr_ord, fill_up)
            if fill_up:
                ptsgrid = grids[j]
            else:
                ptsgrid = np.flipud(grids[j])
            missing = 0
            curr = {"dm41": 0.0, "Ue4sq": 0.0, "Um4sq": 0.0, "Ut4sq": 0.0}
            fill_4 = False
            fill_3 = False
            for currpt in ptsgrid:
                for i in range(3):
                    if fill_up and curr[curr_ord[i]] < currpt[i]:
                        fill_4 = False
                    elif curr[curr_ord[i]] > currpt[i]:
                        fill_3 = False
                for i in range(4):
                    curr[curr_ord[i]] = currpt[i]
                lab = (
                    r"dm41=%.5e " % curr["dm41"]
                    + r"Ue4sq=%.5e " % curr["Ue4sq"]
                    + r"Um4sq=%.5e " % curr["Um4sq"]
                    + r"Ut4sq=%.5e " % curr["Ut4sq"]
                )
                string = "%.5e_%.5e_%.5e_%.5e" % (
                    curr["dm41"],
                    curr["Ue4sq"],
                    curr["Um4sq"],
                    curr["Ut4sq"],
                )
                folder = "grids/%s/OUT/%s/" % (args.gridname, string)
                is_missing = True
                if Neffpoints[string] is not None:
                    is_missing = False
                    if fill_up and Neffpoints[string] > args.fill_up_thres:
                        fill_4 = True
                    elif Neffpoints[string] < args.fill_down_thres:
                        fill_3 = True
                if is_missing:
                    missing += 1
                    if fill_4 or fill_3:
                        if fill_4 and fill_up:
                            print("will fill (4) %s" % folder)
                            if not os.path.exists(folder):
                                os.mkdir(folder)
                            with open("%s/resume.dat" % folder, "w") as _f:
                                _f.write("final w =  1.056\nfinal z =  1.479\n")
                                for i in range(1, 5):
                                    _f.write("dRho_%d  =  0.\n" % i)
                                _f.write("Neff    =  %s\n" % args.fill_up_val)
                            Neffpoints[string] = args.fill_up_val
                            filled += 1
                            missing -= 1
                        if fill_3 and not fill_up:
                            print("will fill (3) %s" % folder)
                            if not os.path.exists(folder):
                                os.mkdir(folder)
                            with open("%s/resume.dat" % folder, "w") as _f:
                                _f.write("final w =  1.097\nfinal z =  1.536\n")
                                for i in range(1, 5):
                                    _f.write("dRho_%d  =  0.\n" % i)
                                _f.write("Neff    =  %s\n" % args.fill_down_val)
                            Neffpoints[string] = args.fill_down_val
                            filled += 1
                            missing -= 1
    print(
        "\n%s: total=%s, missing=%s, filled: %d\n"
        % (args.gridname, len(grid), missing, filled)
    )
    return values, grid, objects


def call_plot(args, gridContent=None):
    def convert_grid(convert_x, convert_y, pts):
        outpts = np.column_stack((pts, np.zeros(np.shape(pts)[0])))
        outpts = np.column_stack((outpts, np.zeros(np.shape(pts)[0])))
        cva = [convert_x, convert_y]
        if "sinsqth14" in cva:
            ix = cva.index("sinsqth14")
            outpts[:, 4 + ix] = pts[:, 1]
        if "sinsqth24" in cva:
            ix = cva.index("sinsqth24")
            outpts[:, 4 + ix] = pts[:, 2] / (1 - pts[:, 1])
        if "sinsqth34" in cva:
            ix = cva.index("sinsqth34")
            outpts[:, 4 + ix] = pts[:, 3] / (1 - pts[:, 1] - pts[:, 2])
        if "sinsq2th_ee" in cva:
            ix = cva.index("sinsq2th_ee")
            outpts[:, 4 + ix] = 4.0 * pts[:, 1] * (1.0 - pts[:, 1])
        if "sinsq2th_mumu" in cva:
            ix = cva.index("sinsq2th_mumu")
            outpts[:, 4 + ix] = 4.0 * pts[:, 2] * (1.0 - pts[:, 2])
        if "sinsq2th_emu" in cva:
            ix = cva.index("sinsq2th_emu")
            outpts[:, 4 + ix] = 4.0 * pts[:, 1] * pts[:, 2]
        return outpts

    # prepare the grid points
    if gridContent:
        mixings, fullgrid, fullobjects = gridContent
    else:
        mixings, fullgrid, fullobjects = call_read(args)
    fullpoints = list(
        map(lambda x: safegetattr(x, "Neff", 0.0) - args.Neff_ref, fullobjects)
    )
    fullpoints = np.asarray(fullpoints)
    cgrid = {}
    if args.par_x == args.par_y:
        print(
            "Cannot plot the same variable in the two axis! %s, %s"
            % (args.par_x, args.par_y)
        )
        return
    fullgrid = convert_grid(args.par_x, args.par_y, fullgrid)
    for a in indexes.keys():
        if getattr(args, "fix_%s" % a) == "False":
            setattr(args, "fix_%s" % a, False)
        if args.par_x == a or args.par_y == a:
            setattr(args, "fix_%s" % a, False)
        if mixings["%s_N" % a] > 1 and getattr(args, "fix_%s" % a) is False:
            cgrid["%s_N" % a] = mixings["%s_N" % a]
            cgrid["%s_min" % a] = mixings["%s_min" % a]
            cgrid["%s_max" % a] = mixings["%s_max" % a]
            cgrid["%s_pts" % a] = mixings["%s_pts" % a]
        else:
            print(
                "fixing %s to %s"
                % (a, mixings["%s_pts" % a][int(getattr(args, "fix_%s" % a))])
            )
            cgrid["%s_N" % a] = 1
            cgrid["%s_min" % a] = mixings["%s_pts" % a][
                int(getattr(args, "fix_%s" % a))
            ]
            cgrid["%s_max" % a] = mixings["%s_pts" % a][
                int(getattr(args, "fix_%s" % a))
            ]
            cgrid["%s_pts" % a] = np.asarray(
                [mixings["%s_pts" % a][int(getattr(args, "fix_%s" % a))]]
            )
    smallgrid = fillGrid(Namespace(**cgrid))
    smallgrid = convert_grid(args.par_x, args.par_y, smallgrid)
    smallpoints = []
    for fgv, pt in zip(fullgrid, fullpoints):
        for ngv in smallgrid:
            if list(fgv) == list(ngv):
                smallpoints.append(pt)
    smallpoints = np.asarray(smallpoints)
    if "sinsq" in args.par_x:
        xv = np.unique(smallgrid[:, 4])
    else:
        xv = mixings["%s_pts" % args.par_x]
    if "sinsq" in args.par_y:
        xv = np.unique(smallgrid[:, 5])
    else:
        yv = mixings["%s_pts" % args.par_y]
    if args.lsn_contours:
        if args.par_y == "dm41" and args.par_x in ["Ue4sq", "sinsqth14", "sinsq2th_ee"]:
            lsn_contours = [
                "/home/gariazzo/data/lsn/globalfit/1801.06469/fig4b-contours/cnt-9973-%d.dat"
                % i
                for i in [1, 2, 3]
            ]
            lsn_contours_convssq2th = args.par_x in ["Ue4sq", "sinsqth14"]
    else:
        lsn_contours = []
        lsn_contours_convssq2th = False
    contourplot(
        xv,
        yv,
        mixings,
        smallpoints,
        fname="grids/%s/plots/%s_%s.pdf" % (args.gridname, args.par_x, args.par_y)
        if args.filename == ""
        else "grids/%s/plots/%s" % (args.gridname, args.filename),
        xlab=labels[args.par_x],
        ylab=labels[args.par_y],
        xlim=[args.xmin, args.xmax] if (args.xmin > 0 and args.xmax > 0) else None,
        ylim=[args.ymin, args.ymax] if (args.ymin > 0 and args.ymax > 0) else None,
        title=args.title,
        bfx=args.bestfit_x,
        bfy=args.bestfit_y,
        bfup=args.bestfit_upper,
        lsn_contours=lsn_contours,
        lsn_contours_convssq2th=lsn_contours_convssq2th,
        colorbar=args.colorbar,
        colorbar_fname=os.path.join(
            "grids", args.gridname, "plots", args.colorbar_fname
        )
        if args.colorbar_fname != ""
        else None,
        textbox=args.textbox if hasattr(args, "textbox") else None,
    )
    print("\nDone!\n\n")
    return


def call_prepare(args):
    if args.gridname in ["plot", "ternary"]:
        print(
            "Sorry, I cannot let you name your grid 'plot' or 'ternary', otherwise I will have problems in the future."
        )
        return
    if not os.path.exists("grids/%s/ini/" % args.gridname):
        os.makedirs("grids/%s/ini/" % args.gridname)
    if not os.path.exists("grids/%s/OUT/" % args.gridname):
        os.makedirs("grids/%s/OUT/" % args.gridname)
    if not os.path.exists("grids/%s/plots/" % args.gridname):
        os.makedirs("grids/%s/plots/" % args.gridname)
    files = list(glob.iglob("grids/%s/ini/*.ini" % args.gridname))
    list(map(lambda x: os.remove(x), files))
    for a in ["dm41", "Ue4sq", "Um4sq", "Ut4sq"]:
        if getattr(args, "%s_N" % a) == 1:
            pmin = getattr(args, "%s_min" % a)
            pmax = getattr(args, "%s_max" % a)
            if pmin == 0.0 or pmax == 0.0:
                setattr(args, "%s_min" % a, 0.5 * (pmin + pmax))
            else:
                setattr(
                    args, "%s_min" % a, 10 ** (0.5 * (np.log10(pmin) + np.log10(pmax)))
                )
            setattr(args, "%s_max" % a, getattr(args, "%s_min" % a))
    if args.ternary:
        for a in ["Um4sq", "Ut4sq"]:
            setattr(args, "%s_min" % a, args.Ue4sq_min)
            setattr(args, "%s_max" % a, args.Ue4sq_max)
            setattr(args, "%s_N" % a, args.Ue4sq_N)
    write_grid_cfg(args)
    grid = fillGrid(args)
    for dm41, Ue4sq, Um4sq, Ut4sq in grid:
        ssq14 = Ue4sq
        ssq24 = Um4sq / (1.0 - Ue4sq)
        ssq34 = Ut4sq / (1.0 - Ue4sq - Um4sq)
        if np.isnan(ssq14):
            ssq14 = 0.0
        if np.isnan(ssq24):
            ssq24 = 0.0
        if np.isnan(ssq34):
            ssq34 = 0.0
        if args.rename:
            oldname = "grids/%s/OUT/%s_%s_%s_%s/" % (
                args.gridname,
                dm41,
                Ue4sq,
                Um4sq,
                Ut4sq,
            )
            newname = "grids/%s/OUT/%.5e_%.5e_%.5e_%.5e/" % (
                args.gridname,
                dm41,
                Ue4sq,
                Um4sq,
                Ut4sq,
            )
            if os.path.exists(oldname):
                os.rename(oldname, newname)
        prep = [
            "grids/%s/ini/%.5e_%.5e_%.5e_%.5e.ini"
            % (args.gridname, dm41, Ue4sq, Um4sq, Ut4sq),
            "grids/%s/OUT/%.5e_%.5e_%.5e_%.5e/"
            % (args.gridname, dm41, Ue4sq, Um4sq, Ut4sq),
            "3+1",
            "damping",
            "--dlsoda_rtol=%s" % args.tolerance,
            "--dlsoda_atol=%s" % args.tolerance,
            "--Nx=%s" % args.Nx,
            "--Ny=%s" % args.Ny,
            "--Nylog=%s" % args.Nylog,
            "--y_cen=%s" % args.y_cen,
            "--x_in=%s" % args.x_in,
            "--default_sterile=None",
            "--dm41=%s" % dm41,
            "--th14=%s" % ssq14,
            "--th24=%s" % ssq24,
            "--th34=%s" % ssq34,
            "--ordering=%s" % args.ordering,
        ]
        if args.no_GL:
            prep.append("--no_GL")
        for s in [
            "save_energy_entropy",
            "save_fd",
            "save_Neff",
            "save_nuDens",
            "save_z",
        ]:
            if getattr(args, s):
                prep.append("--%s" % s)
        parser = prepareIni.setParser()
        rargs = parser.parse_args(prep)
        values = prepareIni.getIniValues(rargs)
        prepareIni.writeIni(rargs.inifile, values)
    print("\nTotal number of points: %s" % len(grid))


def call_read(args):
    write = args.write if hasattr(args, "write") else False
    values, grid = read_grid_cfg(args.gridname)
    objects = []
    missing = 0
    if args.verbose:
        print("\nTotal number of points: %s" % len(grid))
    Neffpoints = {}
    for dm41, Ue4sq, Um4sq, Ut4sq in grid:
        lab = (
            r"dm41=%.5e " % dm41
            + r"Ue4sq=%.5e " % Ue4sq
            + r"Um4sq=%.5e " % Um4sq
            + r"Ut4sq=%.5e " % Ut4sq
        )
        string = "%.5e_%.5e_%.5e_%.5e" % (dm41, Ue4sq, Um4sq, Ut4sq)
        folder = "grids/%s/OUT/%s/" % (args.gridname, string)
        obj = None
        try:
            obj = FortEPiaNORun(
                folder, label=lab, nnu=4, rho=False, verbose=args.verbose
            )
        except (IOError, IndexError):
            if args.verbose:
                print("no %s" % lab)
            missing += 1
        else:
            try:
                Neffpoints[string] = obj.Neff
            except AttributeError:
                Neffpoints[string] = np.nan
                missing += 1
        objects.append(obj)
    if write:
        with open("grids/%s/Neff.dat" % args.gridname, "w") as _f:
            _f.write("#dm41 Ue4sq Um4sq Ut4sq Neff\n")
            for dm41, Ue4sq, Um4sq, Ut4sq in grid:
                s = "%.5e_%.5e_%.5e_%.5e" % (dm41, Ue4sq, Um4sq, Ut4sq)
                _f.write("%s %.5e\n" % (s.replace("_", " "), Neffpoints[s]))
    if args.verbose:
        print("\nMissing or incomplete points: %s\n" % missing)
    else:
        print("\n%s: total=%s, missing=%s\n" % (args.gridname, len(grid), missing))
    return values, grid, objects


def call_run(args):
    print("submitting the grid %s" % args.gridname)
    files = list(glob.iglob("grids/%s/ini/*.ini" % args.gridname))
    if args.failed_only:
        newfiles = []
        for f in files:
            if not os.path.exists(
                f.replace("ini/", "OUT/").replace(".ini", "/resume.dat")
            ):
                if args.remove_existing:
                    shutil.rmtree(
                        f.replace("ini/", "OUT/").replace(".ini", "/"),
                        ignore_errors=True,
                    )
                newfiles.append(f)
        files = newfiles
    current = 0
    files = sortFiles(files, args.from_heaviest)
    if args.verbose:
        print("\n".join(files))
    for i, f in enumerate(files):
        if i >= args.first_index and i < args.last_index:
            if args.local:
                jobcommand = "bin/fortepiano {ini:}".format(ini=f)
            else:
                jobcommand = "clusterlauncher -N {gn:}_{fn:} -n 1 --openmp -q {q:} -w {h:}:{m:}:00 bin/fortepiano {ini:}".format(
                    gn=args.gridname,
                    fn=f.split(os.sep)[-1].replace(".ini", ""),
                    q=args.queue,
                    h=args.walltime_hours,
                    m=args.walltime_minutes,
                    ini=f,
                )
            current += 1
            os.system(jobcommand)
    print("\nTotal number of runs: %s, submitted: %s" % (len(files), current))


def call_ternary(args, gridContent=None):
    if not tern:
        print("ternary not available. Exiting.")
        return
    # prepare the grid points
    if gridContent:
        mixings, fullgrid, fullobjects = gridContent
    else:
        mixings, fullgrid, fullobjects = call_read(args)
    fullpoints = list(
        map(lambda x: safegetattr(x, "Neff", 0.0) - args.Neff_ref, fullobjects)
    )
    fullpoints = np.asarray(fullpoints)
    if (
        mixings["Ue4sq_N"] != mixings["Um4sq_N"]
        or mixings["Ue4sq_N"] != mixings["Ut4sq_N"]
        or mixings["Ue4sq_N"] < 2
    ):
        print(
            "Invalid grid for ternary plot: must have the same number of points (>1)for all the mixing matrix elements"
        )
        return
    cgrid = mixings.copy()
    if args.fix_dm41 >= mixings["dm41_N"]:
        print("dm41 cannot be larger than the number of grid points in dm41!")
        return
    print("fixing dm41 to %s" % (mixings["dm41_pts"][args.fix_dm41]))
    cgrid["dm41_N"] = 1
    cgrid["dm41_min"] = mixings["dm41_pts"][args.fix_dm41]
    cgrid["dm41_max"] = mixings["dm41_pts"][args.fix_dm41]
    cgrid["dm41_pts"] = np.asarray([mixings["dm41_pts"][args.fix_dm41]])
    smallgrid = fillGrid(Namespace(**cgrid))
    smallpoints = []
    for fgv, pt in zip(fullgrid, fullpoints):
        for ngv in smallgrid:
            if list(fgv) == list(ngv):
                smallpoints.append(pt)
    if not mixings["ternary"]:
        smallpoints = np.asarray(smallpoints).reshape(
            (mixings["Ue4sq_N"], mixings["Um4sq_N"], mixings["Ut4sq_N"])
        )

    # prepare the ternary figure
    vmin = 3.0
    vmax = 4.0
    norm = plt.Normalize(vmin=vmin, vmax=vmax)
    scale = mixings["Ue4sq_N"] - (0 if mixings["ternary"] else 1)
    figure, tax = ternary.figure(scale=scale)
    figure.set_size_inches(6, 4)
    tax.gridlines(
        multiple=1,
        linewidth=1,
        linestyle=":",
        horizontal_kwargs={"color": "k", "linestyle": ":"},
        left_kwargs={"color": "k", "linestyle": "--"},
        right_kwargs={"color": "k", "linestyle": "-."},
    )
    tax.boundary(linewidth=2.0)
    colmap = matplotlib.cm.get_cmap("gnuplot")  # tab10, Dark2, Set1, CMRmap

    # call the heatmap
    if not mixings["ternary"]:
        if args.style == "heatmap":

            def function(p):
                """Compute the value of Delta Neff for the simplex points."""
                return smallpoints[
                    int(p[0] * (mixings["Ue4sq_N"] - 1)),
                    int(p[1] * (mixings["Um4sq_N"] - 1)),
                    int(p[2] * (mixings["Ut4sq_N"] - 1)),
                ]

            tax.heatmapf(
                function,
                boundary=True,
                style="t",
                colorbar=False,
                cmap=colmap,
                vmin=vmin,
                vmax=vmax,
            )
        elif args.style == "scatter":
            pts = []
            cs = []
            for pt in ternary.helpers.simplex_iterator(5):
                pts.append(pt)
                cs.append(smallpoints[int(pt[0]), int(pt[1]), int(pt[2])])
            tax.scatter(
                pts,
                s=200,
                c=cs,
                marker="o",
                cmap=colmap,
                norm=norm,
                vmin=vmin,
                vmax=vmax,
                alpha=0.9,
            )
    else:
        if args.style == "heatmap":
            data = dict()
            for n, [i, j, k] in enumerate(
                ternary.helpers.simplex_iterator(mixings["Ue4sq_N"])
            ):
                data[(i, j)] = smallpoints[n]
            tax.heatmap(
                data, style="t", colorbar=False, cmap=colmap, vmin=vmin, vmax=vmax
            )
        elif args.style == "scatter":
            tax.scatter(
                ternary.helpers.simplex_iterator(mixings["Ue4sq_N"]),
                s=200,
                c=smallpoints,
                marker="o",
                cmap=colmap,
                norm=norm,
                vmin=vmin,
                vmax=vmax,
                alpha=0.9,
            )

    # fix tick labels:
    tfontsize = 10
    locations = np.arange(0, scale + 1)
    offset = 0.15
    ticks = [v for v in mixings["Ue4sq_pts"]]
    ax = tax.get_axes()
    for a in ["r", "l", "b"]:
        for index, i in enumerate(locations):
            if a == "r":
                loc1 = (scale - i, i, 0)
                loc2 = (scale - i + offset, i, 0)
                text_location = (scale - i + 1.6 * offset, i + 0.2 - 0.5 * offset, 0)
                tick = ticks[index]
                ha = "left"
                ro = 0
            elif a == "l":
                loc1 = (0, i, 0)
                loc2 = (-offset, i + offset, 0)
                text_location = (-offset, i + 0.3 + 0.2 * offset, 0)
                tick = ticks[-(index + 1)]
                ha = "right"
                ro = 0
            elif a == "b":
                loc1 = (i, 0, 0)
                loc2 = (i, -offset, 0)
                text_location = (i + 0.7 + 0.5 * offset, -0.7 - 0.035 * scale, 0)
                tick = ticks[index]
                ha = "center"
                ro = 0  # 60
            x, y = ternary.project_point(text_location)
            regex = re.match("([0-9.]{3})e([0-9\-]+)", "%.1e" % tick)
            base = float(regex.group(1))
            expon = float(regex.group(2))
            # ticklab = "%.1e"%tick
            ticklab = r"$10^{%.0d}$" % (expon)
            if ("%f" % base).startswith("1.0"):
                ax.text(
                    x,
                    y,
                    ticklab,
                    horizontalalignment=ha,
                    color="k",
                    fontsize=tfontsize,
                    rotation=ro,
                )
                if a == "r":
                    loc1 = (scale - i, i, 0)
                    loc2 = (scale - i + 4 * offset, i, 0)
                elif a == "l":
                    loc1 = (0, i, 0)
                    loc2 = (-4 * offset, i + 4 * offset, 0)
                elif a == "b":
                    loc1 = (i, 0, 0)
                    loc2 = (i, -4 * offset, 0)
            ternary.line(ax, loc1, loc2, color="k")

    # axis labels and colorbar
    tax.clear_matplotlib_ticks()
    tax.get_axes().axis("off")
    fontsize = 14
    offset = 0.15
    ax.text(
        -1,
        scale * 0.75,
        r"$\Delta m^2_{41} = %g$ eV$^2$" % mixings["dm41_pts"][args.fix_dm41],
        # horizontalalignment=ha,
        color="k",
        fontsize=1.2 * fontsize,
    )
    tax.bottom_axis_label(r"$|U_{e4}|^2$", fontsize=fontsize, offset=0.0)
    tax.left_axis_label(r"$|U_{\mu4}|^2$", fontsize=fontsize, offset=offset)
    tax.right_axis_label(r"$|U_{\tau4}|^2$", fontsize=fontsize, offset=offset)
    sm = plt.cm.ScalarMappable(cmap=colmap, norm=norm)
    sm._A = []
    cb = plt.colorbar(sm, extend="both", ticks=[3.1, 3.3, 3.5, 3.7, 3.9])
    cb.set_label(r"$N_{\rm eff}$", fontsize=fontsize)
    plt.tight_layout(rect=(-0.03, -0.035, 1.055, 1.035))

    # save and exit
    tax.savefig(
        "grids/%s/plots/ternary_%s_dm41_%s.pdf"
        % (
            args.gridname,
            "s" if args.style == "scatter" else "h",
            mixings["dm41_pts"][args.fix_dm41],
        )
        if args.filename == ""
        else "grids/%s/plots/%s" % (args.gridname, args.filename)
    )
    print("\nDone!\n\n")
    return


if __name__ == "__main__":
    parser = setParser()
    argsList = sys.argv[1:]
    if (
        argsList.count("plot") > 1
        or argsList.count("ternary") > 1
        or (argsList.count("plot") > 0 and argsList.count("ternary") > 0)
    ):
        if any([x in argsList for x in ["set", "run", "read"]]):
            print(
                "multiple sub-commands can only be used for plotting ('plot' or 'ternary')"
            )
            exit(1)
        if argsList[0] in ["plot", "ternary"]:
            args = parser.parse_args(["read"])
        args = parser.parse_args([argsList[0], "read"])
        gridContent = call_read(args)
        idx = [i for i, x in enumerate(argsList) if x == "plot" or x == "ternary"] + [
            len(argsList)
        ]
        for i in range(len(idx) - 1):
            args = parser.parse_args([argsList[0]] + argsList[idx[i] : idx[i + 1]])
            print(args)
            args.func(args, gridContent)
    else:
        args = parser.parse_args(argsList)
        if args.verbose:
            print(args)
        args.func(args)
