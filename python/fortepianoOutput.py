"""Functions and classes that read the output and help to do plots"""

import os
import re
import traceback

try:
    import matplotlib

    matplotlib.use("agg")
    import matplotlib.pyplot as plt
    from matplotlib.ticker import AutoMinorLocator
except ImportError:
    print("Cannot import matplotlib.pyplot...may raise errors later")
    plt = None
try:
    import numpy as np
except ImportError:
    print("Cannot import numpy...may raise errors later")
    np = None
try:
    from scipy.interpolate import interp1d
    from scipy.integrate import quad
    from scipy.signal import savgol_filter
except ImportError:
    print("Cannot import scipy...may raise errors later")
    interp1d = None
    quad = None

try:
    FileNotFoundError
except NameError:
    FileNotFoundError = IOError

colors = ["r", "g", "b", "k", "c", "m", "#ff9933", "y", "#99ff33"] * 4
styles = ["-", "--", ":", "-."] * 2
markers = [".", "+", "x", "^", "*", "h", "D"]

PISQD15 = np.pi ** 2 / 15.0
ELECTRONMASS_MEV = 0.5109989461


def finalizePlot(
    fname,
    lloc="best",
    title="",
    xlab=None,
    ylab=None,
    xscale=None,
    yscale=None,
    xlim=None,
    ylim=None,
    legcol=1,
    legend=True,
    x_T=False,
    Neff_axes=False,
    tightrect=(-0.025, -0.03, 1.015, 1.03),
):
    """Prepare the final configuration of the plot
    (legends, labels, limits, axes and so on)
    before saving the figure and closing the plot

    Parameters:
        fname: the filename where to save the figure
        lloc (default "best"): the position for the legend
        title (default ""): the title to add to the figure
        xlab, ylab (default None): main labels for the x and y axes
        xscale, yscale (default None): scale ("linear" or "log")
            for the x and y axes
        xlim, ylim (default None):
            the limits to use for the x and y axes
        legcol (default 1): the number of columns to use in the legend
        legend (default True): if True, add a legend to the plot
        x_T (default False): if True, add a twin x axis
            in the upper part of the plot,
            with the conversion betwen x and T
        Neff_axes (default False): if True, the y axis on the left
            will show Neff normalized in the early universe,
            while its twin on the right the final Neff
        tightrect (default (-0.025, -0.03, 1.015, 1.03)):
            the rect parameter to use in tight_plot
    """
    plt.title(title)
    ax = plt.gca()
    if not Neff_axes:
        ax.tick_params(
            "both",
            which="both",
            direction="in",
            left=True,
            right=True,
            top=True,
            bottom=True,
        )
    if legend:
        plt.legend(loc=lloc, ncol=legcol)
    if xlab is not None:
        plt.xlabel(xlab)
    if ylab is not None:
        plt.ylabel(ylab)
    if xscale is not None:
        plt.xscale(xscale)
    if yscale is not None:
        plt.yscale(yscale)
    if xlim is not None:
        plt.xlim(xlim)
    if ylim is not None:
        plt.ylim(ylim)
    if x_T:
        lims = ax.get_xlim()
        ax.set_xscale("log")
        ax.set_xlabel("$x$")
        ax1 = ax.twiny()
        ax1.set_xlim([0.5109989461 / lims[0], 0.5109989461 / lims[1]])
        ax1.set_xscale("log")
        ax1.set_xlabel("$T$ [MeV]")
        ax1.tick_params(
            "both",
            which="both",
            direction="in",
            left=True,
            right=True,
            top=True,
            bottom=False,
        )
    if Neff_axes:
        ax.set_ylabel(r"$N_{\rm eff}^{\rm in}=\frac{8}{7}\frac{\rho_\nu}{\rho_\gamma}$")
        lims = ax.get_ylim()
        ax1 = ax.twinx()
        ax.tick_params(
            "both",
            which="both",
            direction="in",
            left=True,
            right=False,
            labelleft=True,
            labelright=False,
        )
        ax1.tick_params(
            "both",
            which="both",
            direction="in",
            left=False,
            right=True,
            labelleft=False,
            labelright=True,
        )
        ax1.set_ylabel(
            r"$N_{\rm eff}^{\rm now}=\frac{8}{7}\left(\frac{11}{4}\right)^{4/3}\;\frac{\rho_\nu}{\rho_\gamma}$"
        )
        ax1.set_ylim(np.asarray(lims) * (11.0 / 4) ** (4.0 / 3))
        minorLocatorX = AutoMinorLocator(2)
        ax1.yaxis.set_minor_locator(minorLocatorX)
    plt.tight_layout(rect=tightrect)
    try:
        plt.savefig(fname)
    except FileNotFoundError:
        print(traceback.format_exc())
    plt.close()


def stripRepeated(data, ix1, ix2):
    """Strip the repeated points from an output file,
    to avoid steps in the plots

    Parameters:
        data: the np.ndarray with all the plot data
        ix1: the index that corresponds to the x axis in the plot
        ix2: the index that corresponds to the y axis in the plot
    """
    x = []
    y = []
    try:
        for d in data:
            if len(y) == 0:
                x.append(d[ix1])
                y.append(d[ix2])
            if y[-1] != d[ix2]:
                x.append(d[ix1])
                y.append(d[ix2])
        if x[-1] != data[-1][ix1] or y[-1] != data[-1][ix2]:
            x.append(data[-1][ix1])
            y.append(data[-1][ix2])
    except IndexError:
        print(traceback.format_exc())
        return np.asarray([np.nan]), np.asarray([np.nan])
    return np.asarray(x), np.asarray(y)


class FortEPiaNORun:
    """Class that reads the output of FortEPiaNO and helps to do plots,
    compute integrals of the density matrix or other things
    """

    def __init__(
        self,
        folder,
        nnu=3,
        full=True,
        label="",
        plots=False,
        verbose=True,
        deltas=False,
    ):
        """Read the entire output of FortEPiaNO from a specific folder.
        It will ignore non-existing files and store all the available
        information for further processing (plots, ...)

        Parameters:
            folder: the name of the folder to consider
            nnu (default 3): number of neutrinos to consider
                (if more rows/columns of the density matrix exist,
                they will be ignored)
            full (default True): if True, read also all the off-diagonal
                density matrix elements, otherwise ignore them
                (to save time if not needed in the plots, for example)
            label (default ""): a label to assign
                to the current FortEPiaNO point in the plots
            plots (default False): if True, produce a series of plots
                after having read all the files
            verbose (default True): if True, print more error messages
                (e.g. when the folder is not found or w is not saved)
            deltas (default False): if True, print the relative variation
                of energy and number density for each neutrino
        """
        self.folder = folder
        self.full = full
        self.label = label
        self.verbose = verbose
        self.nnu = nnu
        self.hasBBN = False
        self.lowReheating = False
        if not os.path.exists(folder):
            if verbose:
                print("non-existing folder: %s" % folder)
            return
        self.zCol = 1
        self.Trhini = None
        try:
            with open("%s/ini.log" % folder) as _ini:
                self.ini = _ini.read()
        except FileNotFoundError:
            self.ini = ""
        else:
            self.ini = self.ini.replace("\n", " ")
            try:
                self.Trhini = float(
                    re.search("Trh =[ ]*([E+\-\d.]*)", self.ini).group(1)
                )
            except (AttributeError, ValueError):
                pass
            else:
                if self.Trhini > 0:
                    self.lowReheating = True
                    self.zCol = 2
        try:
            fdy = np.loadtxt("%s/fd.dat" % folder)
        except (IOError, OSError):
            self.yv = np.nan
            self.fd = np.nan
        else:
            self.yv = fdy[:, 0]
            self.fd = fdy[:, 1]
        try:
            self.zdat = np.loadtxt("%s/z.dat" % folder)
        except (IOError, OSError):
            self.zdat = np.asarray([[np.nan, np.nan, np.nan]])
        try:
            self.Neffdat = np.loadtxt("%s/Neff.dat" % folder)
        except (IOError, OSError):
            self.Neffdat = np.asarray(
                [[np.nan, np.nan, np.nan, np.nan]]
                if self.lowReheating
                else [[np.nan, np.nan, np.nan]]
            )
        try:
            self.endens = np.loadtxt("%s/energyDensity.dat" % folder)
        except (IOError, OSError):
            self.endens = np.nan
        else:
            try:
                self.endens[:, :]
            except IndexError:
                pass
            else:
                if deltas:
                    self.delta_ed = [
                        (
                            self.endens[-1, self.zCol + 4 + i]
                            - self.endens[0, self.zCol + 4 + i]
                        )
                        / self.endens[0, self.zCol + 4 + i]
                        * 100
                        for i in range(nnu)
                    ]
                    self.tot_delta_ed = (
                        (
                            np.sum(self.endens[-1, self.zCol + 4 : self.zCol + 4 + nnu])
                            - np.sum(
                                self.endens[0, self.zCol + 4 : self.zCol + 4 + nnu]
                            )
                        )
                        / np.sum(self.endens[0, self.zCol + 4 : self.zCol + 4 + nnu])
                        * 100
                    )
                    print(
                        "delta energy density:\t"
                        + "\t".join(
                            [
                                "nu%d: %f%%"
                                % (
                                    i + 1,
                                    self.delta_ed[i],
                                )
                                for i in range(nnu)
                            ]
                        )
                    )
        try:
            self.entropy = np.loadtxt("%s/entropy.dat" % folder)
        except (IOError, OSError):
            self.entropy = np.nan
        try:
            self.number = np.loadtxt("%s/numberDensity.dat" % folder)
        except (IOError, OSError):
            self.number = np.nan
        else:
            try:
                self.number[:, :]
            except IndexError:
                pass
            else:
                if deltas:
                    self.delta_nd = [
                        (
                            self.number[-1, self.zCol + 4 + i]
                            - self.number[0, self.zCol + 4 + i]
                        )
                        / self.number[0, self.zCol + 4 + i]
                        * 100
                        for i in range(nnu)
                    ]
                    self.tot_delta_nd = (
                        (
                            np.sum(self.number[-1, self.zCol + 4 : self.zCol + 4 + nnu])
                            - np.sum(
                                self.number[0, self.zCol + 4 : self.zCol + 4 + nnu]
                            )
                        )
                        / np.sum(self.number[0, self.zCol + 4 : self.zCol + 4 + nnu])
                        * 100
                    )
                    print(
                        "delta number density:\t"
                        + "\t".join(
                            [
                                "nu%d: %f%%" % (i + 1, self.delta_nd[i])
                                for i in range(nnu)
                            ]
                        )
                    )
        self.rho = np.asarray([[[None, None] for i in range(nnu)] for j in range(nnu)])
        self.rhoM = np.asarray([[[None, None] for i in range(nnu)] for j in range(nnu)])
        try:
            with open("%s/resume.dat" % folder) as _f:
                self.resume = _f.read()
        except FileNotFoundError:
            self.resume = ""
            self.hasResume = False
        else:
            self.resume = self.resume.replace("\n", " ")
            self.hasResume = True
        if self.hasResume:
            try:
                self.Neff = float(re.search("Neff[ =]*([-\d.]*)", self.resume).group(1))
            except (AttributeError, ValueError):
                self.Neff = np.nan
            try:
                self.wfin = float(
                    re.search("final w[ =]*([-\d.]*)", self.resume).group(1)
                )
            except (AttributeError, ValueError):
                if verbose:
                    print("cannot read w from resume.dat")
                self.wfin = np.nan
            try:
                self.zfin = float(
                    re.search("final z[ =]*([-\d.]*)", self.resume).group(1)
                )
            except (AttributeError, ValueError):
                self.zfin = np.nan
            try:
                self.Trh = float(re.search("Trh[ =]*([-\d.]*)", self.resume).group(1))
            except (AttributeError, ValueError):
                self.Trh = np.nan
            else:
                if not np.isclose(float(self.Trh), float(self.Trhini)):
                    raise ValueError(
                        "Trh from ini.log (%e) and from resume.dat (%e) differ."
                        % (self.Trh, self.Trhini)
                    )
        self.deltarhofin = []
        for i in range(self.nnu):
            if self.hasResume:
                try:
                    self.deltarhofin.append(
                        float(
                            re.search(
                                "dRho_%s[ =]*([-\d.]*)" % (i + 1),
                                self.resume,
                            ).group(1)
                        )
                    )
                except (AttributeError, IndexError):
                    self.deltarhofin.append(np.nan)
            try:
                self.rho[i, i, 0] = np.loadtxt("%s/nuDens_diag%d.dat" % (folder, i + 1))
            except (IOError, OSError):
                self.rho[i, i, 0] = np.nan
            if full:
                for j in range(i + 1, self.nnu):
                    try:
                        self.rho[i, j, 0] = np.loadtxt(
                            "%s/nuDens_nd_%d%d_re.dat" % (folder, i + 1, j + 1)
                        )
                    except (IOError, OSError):
                        self.rho[i, j, 0] = np.nan
                    try:
                        self.rho[i, j, 1] = np.loadtxt(
                            "%s/nuDens_nd_%d%d_im.dat" % (folder, i + 1, j + 1)
                        )
                    except (IOError, OSError):
                        self.rho[i, j, 1] = np.nan
            try:
                self.rhoM[i, i, 0] = np.loadtxt(
                    "%s/nuDens_mass%d.dat" % (folder, i + 1)
                )
            except (IOError, OSError):
                self.rhoM[i, i, 0] = np.nan
        self.prepareBBN()
        self.printTableLine()
        if plots:
            self.doAllPlots()

    def prepareBBN(self):
        """Read BBN files and prepare output for PArthENoPE"""
        folder = self.folder
        if (
            os.path.exists("%s/BBN.dat" % folder)
            and os.path.exists("%s/rho_tot.dat" % folder)
            and os.path.exists("%s/nuDens_diag1_BBN.dat" % folder)
            and not os.path.exists("%s/parthenope.dat" % folder)
            and not os.path.exists("%s/parthenope_yi.dat" % folder)
            and not os.path.exists("%s/parthenope_rhoee.dat" % folder)
        ):
            self.hasBBN = True
            # read main quantities and prepare N function
            self.bbn = np.loadtxt("%s/BBN.dat" % folder)
            try:
                self.bbn[:, 3]
                xvec = self.bbn[:, 0]
                assert np.allclose(xvec, self.zdat[:, 0])
                assert np.allclose(
                    self.bbn[:, 1], self.zdat[:, 2 if self.lowReheating else 1]
                )
            except (TypeError, IndexError):
                print("Error in the structure of the BBN.dat file: cannot find columns")
                self.hasBBN = False
            else:
                # read x values at which (rho_rad >= 0.99 rho_tot) for the first time
                self.summedrhos = np.loadtxt("%s/rho_tot.dat" % folder)
                x99 = 0.0
                try:
                    rho_tot = self.summedrhos[:, 2 if self.lowReheating else 1]
                    rho_rad = self.summedrhos[:, 1]
                except (TypeError, IndexError):
                    print("Error in the structure of rho_tot.dat: cannot find columns")
                    self.hasBBN = False
                else:
                    x99 = next(
                        x for x, t, r in zip(xvec, rho_tot, rho_rad) if r >= 0.99 * t
                    )
                    self.filter99 = np.where(xvec >= x99)
                    self.parthenope = np.c_[
                        self.bbn[self.filter99],
                        np.clip(self.drhonu_dx_savgol[self.filter99], 1e-16, None),
                        rho_rad[self.filter99],
                        rho_tot[self.filter99],
                    ]
                    self.parthenope_cols = [
                        "x",
                        "z",
                        "w",
                        "rhobarnu",
                        "drhobarnu_dx",
                        "rho_rad",
                        "rho_tot",
                    ]
                    try:
                        data = np.loadtxt("%s/nuDens_diag1_BBN.dat" % folder)
                        data[:, 2:][self.filter99]
                    except (IndexError, TypeError):
                        print(
                            "Error in the structure of nuDens_diag1_BBN.dat: cannot find columns"
                        )
                    else:
                        try:
                            np.savetxt(
                                "%s/parthenope.dat" % folder,
                                self.parthenope,
                                fmt="%15.7e",
                                header="".join(
                                    ["%16s" % s for s in self.parthenope_cols]
                                )[3:],
                            )
                            np.savetxt(
                                "%s/parthenope_yi.dat" % folder, self.yv, fmt="%15.7e"
                            )
                            np.savetxt(
                                "%s/parthenope_rhoee.dat" % folder,
                                data[:, 2:][self.filter99],
                                fmt="%15.7e",
                            )
                        except IOError:
                            print("Cannot write the outputs for PArtENoPE!")

    @property
    def x(self):
        """return x"""
        return self.zdat[:, 0]

    @property
    def z(self):
        """return z"""
        return self.zdat[:, 2 if self.lowReheating else 1]

    @property
    def Tgamma(self):
        """return T_gamma"""
        return self.z * ELECTRONMASS_MEV / self.x

    @property
    def t(self):
        """return t, if available"""
        return self.zdat[:, 1] if self.lowReheating else np.nan

    @property
    def rhonu(self):
        """return rhonu (total neutrino energy density)"""
        return self.bbn[:, 3] if self.hasBBN else np.nan

    @property
    def drhonu_dx(self):
        """return drhonu/dx"""
        return np.gradient(self.rhonu, self.x) if self.hasBBN else np.nan

    @property
    def drhonu_dx_savgol(self):
        """return a filtered version of drhonu/dx"""
        return (
            savgol_filter(np.clip(self.drhonu_dx, 1e-10, None), 51, 1)
            if self.hasBBN
            else np.nan
        )

    @property
    def N_func(self):
        """return the N function used in PArthENoPE"""
        return self.x * self.drhonu_dx / self.z ** 4 if self.hasBBN else np.nan

    @property
    def N_savgol(self):
        """return a filtered version of the N function used in PArthENoPE"""
        return (
            savgol_filter(np.clip(self.N_func, 1e-11, None), 75, 1)
            if self.hasBBN
            else np.nan
        )

    def interpolateRhoIJ(self, i1, i2, y, ri=0, y2=False, mass=False):
        """Interpolate any entry of the density matrix at a given y,
        and return its value for all the saved x points.
        Repeated points with the same f(x, y) at different x
        will be reported only once, at the earlier x.

        Parameters:
            i1: row index of the density matrix to consider
            i2: column index of the density matrix to consider
            y: the y value at which to interpolate
            ri (default 0): it 0, use real part, if 1 the imaginary one
                (only for off-diagonal entries of the density matrix)
            y2 (default False): if True, multiply the output by y**2
            mass (default False): if True, use the density matrix
                in the mass basis

        Output:
            two 1D np.ndarrays of the same size,
            one with the x points and
            one with the values of the function at each y
        """
        if mass:
            rho = self.rhoM
        else:
            rho = self.rho
        xv = []
        yv = []
        prevy = 0
        for i, x in enumerate(rho[i1, i2, ri][:, 0]):
            fy = interp1d(self.yv, rho[i1, i2, ri][i, 1:])
            cy = fy(y) * (y ** 2 if y2 else 1.0)
            if cy != prevy:
                prevy = cy
                yv.append(prevy)
                xv.append(x)
        if x != xv[-1] or cy != yv[-1]:
            xv.append(x)
            yv.append(cy)
        return np.asarray(xv), np.asarray(yv)

    def interpolateRhoIJ_x(self, i1, i2, x, ri=0, y2=False, mass=False):
        """Interpolate any entry of the density matrix at a given x,
        and return its value for all the y grid points

        Parameters:
            i1: row index of the density matrix to consider
            i2: column index of the density matrix to consider
            x: the x value at which to interpolate
            ri (default 0): it 0, use real part, if 1 the imaginary one
                (only for off-diagonal entries of the density matrix)
            y2 (default False): if True, multiply the output by y**2
            mass (default False): if True, use the density matrix
                in the mass basis

        Output:
            two 1D np.ndarrays of the same size,
            one with the y grid and
            one with the values of the function at each y
        """
        if mass:
            rho = self.rhoM
        else:
            rho = self.rho
        ov = []
        for i, y in enumerate(self.yv):
            fx = interp1d(
                rho[i1, i2, ri][:, 0],
                rho[i1, i2, ri][:, i + 1] * (y ** 2 if y2 else 1.0),
            )
            ov.append(fx(x))
        return self.yv, np.asarray(ov)

    def printTableLine(self):
        """Print a summary of the current point in a latex table format.
        It will contain the label, the final z, all the deltaRho
        computed for the diagonal elements, and the final Neff.
        If the run was not complete, print just the label and a message
        on the last step that was saved
        """
        try:
            self.hasResume
            self.label, self.Neff, self.zfin, self.zdat[-1]
            [self.deltarhofin[i] for i in range(self.nnu)]
        except (AttributeError, IndexError):
            print(traceback.format_exc())
            return
        if not self.verbose:
            return
        if self.hasResume:
            deltastr = ""
            for i in range(self.nnu):
                deltastr += "{:.5f} & ".format(self.deltarhofin[i])
            print(
                "{lab:<15s} & {zfin:.5f} & {deltastr:s}{Neff:.5f}\\\\".format(
                    lab=self.label, Neff=self.Neff, zfin=self.zfin, deltastr=deltastr
                )
            )
        else:
            print(
                "{lab:<35s} \t\t\t\tnot finished, currently on x={x:}".format(
                    lab=self.label, x=self.zdat[-1]
                )
            )

    def plotFD(self, ls="-", lc="k", lab=None, rescale=1.0, fac=1.0):
        """Plot a Fermi-Dirac distribution in the adopted momentum grid.
        It may be rescaled appropriately

        Parameters:
            ls (default "-"): the line style
            lc (default "k"): the line color
            lab (default None): if not None, the line label
            rescale (default 1.0): a rescaling factor of the temperature
                for the Fermi-Dirac distribution
                {enters as f in 1/[exp(y/f) + 1] }
            fac (default 1.0): a constant normalization factor
                for the function that is plotted
        """
        try:
            self.fd, self.yv
        except AttributeError:
            print(traceback.format_exc())
            return
        if rescale != 1.0:
            fd = self.fd * (np.exp(self.yv) + 1.0) / (np.exp(self.yv / rescale) + 1.0)
        else:
            fd = self.fd
        plt.plot(
            self.yv,
            fac * fd,
            label=self.label if lab is None else lab,
            ls=ls,
            marker=".",
            c=lc,
        )
        plt.xscale("log")
        plt.yscale("log")
        plt.xlabel("$y$")
        plt.ylabel(r"$y^2 f(y)$")

    def plotZ(self, ls="-", lc="k", lab=None):
        """Plot z as a function of x

        Parameters:
            ls (default "-"): the line style
            lc (default "k"): the line color
            lab (default None): if not None, the line label
        """
        try:
            self.zdat
        except AttributeError:
            print(traceback.format_exc())
            return
        plt.plot(
            *stripRepeated(self.zdat, 0, 1),
            label=self.label if lab is None else lab,
            ls=ls,
            c=lc
        )
        plt.xscale("log")
        plt.xlabel("$x$")
        plt.ylabel(r"$z$")

    def plotW(self, ls="-", lc="k", lab=None):
        """Plot w as a function of x

        Parameters:
            ls (default "-"): the line style
            lc (default "k"): the line color
            lab (default None): if not None, the line label
        """
        try:
            self.zdat[0, 2]
        except (AttributeError, IndexError):
            print("w is not in z.dat")
            return
        plt.plot(
            *stripRepeated(self.zdat, 0, 2),
            label=self.label if lab is None else lab,
            ls=ls,
            c=lc
        )
        plt.xscale("log")
        plt.xlabel("$x$")
        plt.ylabel(r"$w$")

    def plotZoverW(self, ls="-", lc="k", lab=None):
        """Plot the ratio w/z as a function of x

        Parameters:
            ls (default "-"): the line style
            lc (default "k"): the line color
            lab (default None): if not None, the line label
        """
        try:
            self.zdat[0, 2]
        except (AttributeError, IndexError):
            print("w is not in z.dat")
            return
        plt.plot(
            *stripRepeated(np.asarray([[x[0], x[1] / x[2]] for x in self.zdat]), 0, 1),
            label=self.label if lab is None else lab,
            ls=ls,
            c=lc
        )
        plt.xscale("log")
        plt.xlabel("$x$")
        plt.ylabel(r"$z/w$")

    def plotDeltaZ(self, ref, ls="-", lc="k", lab=None):
        """Plot the difference in z with respect
        to a different FortEPiaNO point, as a function of x

        Parameters:
            ref: a different FortEPiaNORun instance
            ls (default "-"): the line style
            lc (default "k"): the line color
            lab (default None): if not None, the line label
        """
        try:
            self.zdat, ref.zdat
        except AttributeError:
            print(traceback.format_exc())
            return
        mex, mey = stripRepeated(self.zdat, 0, 1)
        mef = interp1d(mex, mey)
        refx, refy = stripRepeated(ref.zdat, 0, 1)
        reff = interp1d(refx, refy)
        plt.plot(
            mex,
            reff(mex) - mef(mex),
            label=self.label if lab is None else lab,
            ls=ls,
            c=lc,
        )
        plt.xscale("log")
        plt.xlabel("$x$")
        plt.ylabel(r"$z-z_{\rm ref}$")

    def plotRhoDiag(self, inu, iy, ls, lc="k", mass=False):
        """Plot one diagonal element of the density matrix
        at a given point (index iy) in the momentum grid,
        as a function of x

        Parameters:
            inu: diagonal index of the density matrix entry to consider
            iy: the index of the requested momentum in the momentum grid
            ls: the line style
            lc (default "k"): the line color
            mass (default False): if True, use the density matrix
                in the mass basis
        """
        try:
            if mass:
                rho = self.rhoM
            else:
                rho = self.rho
        except AttributeError:
            print(traceback.format_exc())
            return
        try:
            rho[inu, inu, 0][:]
        except TypeError:
            print(traceback.format_exc())
            return
        plt.plot(
            *stripRepeated(rho[inu, inu, 0], 0, iy),
            label=r"%s $\alpha$=%d" % (self.label, inu + 1),
            ls=ls,
            c=lc
        )
        plt.xscale("log")
        plt.xlabel("$x$")
        plt.ylabel(r"$\rho_{\alpha\alpha}$")

    def plotdRhoDiag(self, inu, iy, ls, lc="k", mass=False):
        """Plot the x derivative (np.gradient)
        of one diagonal element of the density matrix
        at a given point (index iy) in the momentum grid,
        as a function of x

        Parameters:
            inu: diagonal index of the density matrix entry to consider
            iy: the index of the requested momentum in the momentum grid
            ls: the line style
            lc (default "k"): the line color
            mass (default False): if True, use the density matrix
                in the mass basis
        """
        try:
            if mass:
                rho = self.rhoM
            else:
                rho = self.rho
        except AttributeError:
            print(traceback.format_exc())
            return
        try:
            rho[inu, inu, 0][:]
        except TypeError:
            print(traceback.format_exc())
            return
        dijrex, dijrey = stripRepeated(rho[inu, inu, 0], 0, iy)
        plt.plot(
            dijrex,
            np.gradient(dijrey, dijrex),
            ls=ls,
            c=lc,
            label=r"%s $\alpha$=%d" % (self.label, inu + 1),
        )
        plt.xscale("log")
        plt.xlabel("$x$")
        plt.ylabel(r"$d\rho_{\alpha\alpha}/dx$")

    def plotRhoOffDiag(self, i1, i2, iy, lc="k", im=True, mass=False):
        """Plot one off-diagonal element of the density matrix
        at a given point (index iy) in the momentum grid,
        as a function of x

        Parameters:
            i1: row index of the density matrix to consider
            i2: column index of the density matrix to consider
            iy: the index of the requested momentum in the momentum grid
            lc (default "k"): the line color
            im (default True): add plot of the imaginary part
            mass (default False): if True, use the density matrix
                in the mass basis
        """
        if not hasattr(self, "full") or not self.full:
            print("no offdiagonal loaded")
            return
        try:
            if mass:
                rho = self.rhoM
            else:
                rho = self.rho
        except AttributeError:
            print(traceback.format_exc())
            return
        try:
            rho[i1, i2, 0][:]
        except TypeError:
            print(traceback.format_exc())
            return
        plt.plot(
            *stripRepeated(rho[i1, i2, 0], 0, iy),
            ls="-",
            c=lc,
            label=r"%s $\alpha\beta$=%d%d re" % (self.label, i1 + 1, i2 + 1)
        )
        if im:
            plt.plot(
                *stripRepeated(rho[i1, i2, 1], 0, iy),
                ls=":",
                c=lc,
                label=r"%s $\alpha\beta$=%d%d im" % (self.label, i1 + 1, i2 + 1)
            )
        plt.xscale("log")
        plt.xlabel("$x$")
        plt.ylabel(r"$\rho_{\alpha\beta}$")

    def plotdRhoOffDiag(self, i1, i2, iy, lc="k", im=True, mass=False):
        """Plot the x derivative (np.gradient)
        of one off-diagonal element of the density matrix
        at a given point (index iy) in the momentum grid,
        as a function of x

        Parameters:
            i1: row index of the density matrix to consider
            i2: column index of the density matrix to consider
            iy: the index of the requested momentum in the momentum grid
            lc (default "k"): the line color
            im (default True): add plot of the imaginary part
            mass (default False): if True, use the density matrix
                in the mass basis
        """
        if not hasattr(self, "full") or not self.full:
            print("no offdiagonal loaded")
            return
        try:
            if mass:
                rho = self.rhoM
            else:
                rho = self.rho
            dijrex, dijrey = stripRepeated(rho[i1, i2, 0], 0, iy)
        except (AttributeError, TypeError):
            print(traceback.format_exc())
            return
        plt.plot(
            dijrex,
            np.gradient(dijrey, dijrex),
            ls="-",
            c=lc,
            label=r"%s $\alpha\beta$=%d%d re" % (self.label, i1 + 1, i2 + 1),
        )
        if im:
            dijimx, dijimy = stripRepeated(rho[i1, i2, 1], 0, iy)
            plt.plot(
                dijimx,
                np.gradient(dijimy, dijimx),
                ls=":",
                c=lc,
                label=r"%s $\alpha\beta$=%d%d im" % (self.label, i1 + 1, i2 + 1),
            )
        plt.xscale("log")
        plt.xlabel("$x$")
        plt.ylabel(r"$d\rho_{\alpha\beta}/dx$")

    def plotRhoFin(
        self, i1, i2=None, ri=0, ls="-", lc="k", y2=False, lab=None, mass=False
    ):
        """Plot the y dependence of an element of the density matrix
        at the final x

        Parameters:
            i1: diagonal element of the density matrix to consider,
                or row index of the density matrix element
            i2 (default None): if not None, the column index
                of the requested density matrix element.
                If None, fix i2 = i1 and use a diagonal element
            ri (default 0): it 0, use real part, if 1 the imaginary one
                (only for off-diagonal entries of the density matrix)
            ls (default "-"): the line style
            lc (default "k"): the line color
            y2 (default False): if True,
                multiply the diagonal elements times y**2
            lab (default None): if not None, the line label
            mass (default False): if True, use the density matrix
                in the mass basis
        """
        if i2 is None:
            i2 = i1
        if ri not in [0, 1]:
            ri = 0
        try:
            if mass:
                rho = self.rhoM
            else:
                rho = self.rho
        except AttributeError:
            print(traceback.format_exc())
            return
        try:
            rho[i1, i2, ri][-1, 1:]
        except TypeError:
            print(traceback.format_exc())
            return
        label = (
            r"%s $\alpha\beta$=%d%d %s"
            % (self.label, i1 + 1, i2 + 1, "re" if ri == 0 else "im")
            if lab is None
            else lab
        )
        fyv = self.yv ** 2 * rho[i1, i2, ri][-1, 1:] if y2 else rho[i1, i2, ri][-1, 1:]
        plt.plot(self.yv, fyv, ls=ls, c=lc, label=label)
        plt.xlabel("$y$")
        plt.ylabel(r"$%s\rho_{\alpha\beta}^{\rm fin}(y)$" % ("y^2" if y2 else ""))

    def plotRhoX(
        self,
        i1,
        x,
        i2=None,
        ri=0,
        ls="-",
        lc="k",
        lab="",
        y2=False,
        mass=False,
        divide_by=1.0,
        divide_fd=False,
    ):
        """Plot the y dependence of an element of the density matrix
        at a given x

        Parameters:
            i1: row index of the density matrix entry to consider
            x: the x value at which the density matrix is interpolated
            i2 (default None): if not None, the column index
                of the requested density matrix element.
                If None, fix i2 = i1 and use a diagonal element
            ri (default 0): it 0, use real part, if 1 the imaginary one
                (only for off-diagonal entries of the density matrix)
            ls (default "-"): the line style
            lc (default "k"): the line color
            lab (default ""): if not empty, it will be used as line label
            y2 (default False): if True,
                multiply the diagonal elements times y**2
            mass (default False): if True, use the density matrix
                in the mass basis
            divide_by (default 1.0): divide the rho values by the given
                float or array
            divide_fd (default False): if True,
                divide by self.fd
        """
        if i2 is None:
            i2 = i1
        if ri not in [0, 1]:
            ri = 0
        try:
            interp = self.interpolateRhoIJ_x(i1, i2, x, ri, y2=y2, mass=mass)
        except (AttributeError, TypeError):
            print(traceback.format_exc())
            return
        xv, yv = interp
        plt.plot(
            xv,
            yv / divide_by / (self.fd if divide_fd else 1.0),
            ls=ls,
            c=lc,
            label=r"%s $\alpha\beta$=%d%d %s x=%f"
            % (self.label, i1 + 1, i2 + 1, "re" if ri == 0 else "im", x)
            if not lab
            else lab,
        )
        plt.xlabel("$y$")
        plt.ylabel(r"$%s\rho_{\alpha\beta}(y)$" % ("y^2" if y2 else ""))

    def plotRhoDiagY(
        self, inu, y, ls, lc="k", lab=None, y2=False, mass=False, divide_by=1.0
    ):
        """Plot one diagonal element of the density matrix at a given y
        as a function of x

        Parameters:
            inu: diagonal element of the density matrix to consider
            y: the y value at which the density matrix is interpolated
            ls: the line style
            lc (default "k"): the line color
            lab (default None): if not None, the line label
            y2 (default False): if True,
                multiply the diagonal elements times y**2
            mass (default False): if True, use the density matrix
                in the mass basis
            divide_by (default 1.0): divide the rho values by the given
                float or array
        """
        try:
            x, yv = self.interpolateRhoIJ(inu, inu, y, ri=0, mass=mass)
        except (AttributeError, TypeError):
            print(traceback.format_exc())
            return
        label = lab if lab is not None else r"%s $\alpha$=%d" % (self.label, inu + 1)
        plt.plot(
            x,
            np.asarray(yv) * (y ** 2 if y2 else 1.0) / divide_by,
            label=label,
            ls=ls,
            c=lc,
        )
        plt.xscale("log")
        plt.xlabel("$x$")
        plt.ylabel(r"$%s\rho_{\alpha\alpha}$" % ("y^2" if y2 else ""))

    def plotdRhoDiagY(self, inu, y, ls, lc="k", lab=None, y2=False, mass=False):
        """Plot the x derivative (np.gradient)
        of one diagonal element of the density matrix at a given y
        as a function of x

        Parameters:
            inu: diagonal element of the density matrix to consider
            y: the y value at which the density matrix is interpolated
            ls: the line style
            lc (default "k"): the line color
            lab (default None): if not None, the line label
            y2 (default False): if True,
                multiply the diagonal elements times y**2
            mass (default False): if True, use the density matrix
                in the mass basis
        """
        try:
            x, yv = self.interpolateRhoIJ(inu, inu, y, ri=0, mass=mass)
        except (AttributeError, TypeError):
            print(traceback.format_exc())
            return
        label = lab if lab is not None else r"%s $\alpha$=%d" % (self.label, inu + 1)
        plt.plot(
            x,
            np.gradient(np.asarray(yv) * (y ** 2 if y2 else 1.0), x),
            label=label,
            ls=ls,
            c=lc,
        )
        plt.xscale("log")
        plt.xlabel("$x$")
        plt.ylabel(r"$d%s\rho_{\alpha\alpha}/dx$" % ("y^2" if y2 else ""))

    def plotRhoOffDiagY(self, i1, i2, y, lc="k", ls="-", im=True, lab=None, mass=False):
        """Plot one off-diagonal element of the density matrix
        at a given y as a function of x

        Parameters:
            i1: row index of the density matrix to consider
            i2: column index of the density matrix to consider
            y: the y value at which the density matrix is interpolated
            lc (default "k"): the line color
            ls (default "-"): the line style
            im (default True): add plot of the imaginary part
            lab (default None): if not None, the line label
            mass (default False): if True, use the density matrix
                in the mass basis
        """
        if not hasattr(self, "full") or not self.full:
            print("no offdiagonal loaded")
            return
        try:
            interp = self.interpolateRhoIJ(i1, i2, y, ri=0, mass=mass)
        except (AttributeError, TypeError):
            print(traceback.format_exc())
            return
        plt.plot(
            *interp,
            ls=ls,
            c=lc,
            label=r"%s $\alpha\beta$=%d%d re" % (self.label, i1 + 1, i2 + 1)
            if lab is None
            else lab
        )
        if im:
            try:
                interpIm = self.interpolateRhoIJ(i1, i2, y, ri=1)
            except AttributeError:
                print(traceback.format_exc())
                return
            plt.plot(
                *interpIm,
                ls=":",
                c=lc,
                label=r"%s $\alpha\beta$=%d%d im" % (self.label, i1 + 1, i2 + 1)
                if lab is None
                else lab
            )
        plt.xscale("log")
        plt.xlabel("$x$")
        plt.ylabel(r"$\rho_{\alpha\beta}$")

    def plotdRhoOffDiagY(
        self, i1, i2, y, lc="k", ls="-", im=True, lab=None, mass=False
    ):
        """Plot the x derivative (np.gradient)
        of one off-diagonal element of the density matrix at a given y
        as a function of x

        Parameters:
            i1: row index of the density matrix to consider
            i2: column index of the density matrix to consider
            y: the y value at which the density matrix is interpolated
            lc (default "k"): the line color
            ls (default "-"): the line style
            im (default True): add plot of the imaginary part
            lab (default None): if not None, the line label
            mass (default False): if True, use the density matrix
                in the mass basis
        """
        if not hasattr(self, "full") or not self.full:
            print("no offdiagonal loaded")
            return
        try:
            dijrex, dijrey = self.interpolateRhoIJ(i1, i2, y, ri=0, mass=mass)
        except (AttributeError, TypeError):
            print(traceback.format_exc())
            return
        try:
            plt.plot(
                dijrex,
                np.gradient(dijrey, dijrex),
                ls=ls,
                c=lc,
                label=r"%s $\alpha\beta$=%d%d re" % (self.label, i1 + 1, i2 + 1)
                if lab is None
                else lab,
            )
        except IndexError:
            pass
        if im:
            try:
                dijimx, dijimy = self.interpolateRhoIJ(i1, i2, y, ri=1, mass=mass)
            except AttributeError:
                print(traceback.format_exc())
                return
            try:
                plt.plot(
                    dijimx,
                    np.gradient(dijimy, dijimx),
                    ls=":",
                    c=lc,
                    label=r"%s $\alpha\beta$=%d%d im" % (self.label, i1 + 1, i2 + 1)
                    if lab is None
                    else lab,
                )
            except IndexError:
                pass
        plt.xscale("log")
        plt.xlabel("$x$")
        plt.ylabel(r"$d\rho_{\alpha\beta}/dx$")

    def plotNeff(
        self,
        lc="k",
        ls="-",
        lab=None,
        nefflims=[0.5, 4.5],
        axes=True,
        endensx=0,
        endensgamma=2,
        endensnu0=5,
        useT=False,
    ):
        """Plot the evolution of Neff as a function of x.
        If not available from Neff.dat,
        compute it directly from the density matrix

        Parameters:
            lc (default "k"): the line color
            ls (default "-"): the line style
            lab (default None): if not None, the line label
            nefflims (default [0.5, 4.5]): range for the y axis
            axes (default True): if True, create a twin y axis
                with Neff^now in addition to the one with Neff^initial
            endensx (default 0): index of
                the x column in self.endens
            endensgamma (default 2): index of
                the photon energy density column in self.endens
            endensnu0 (default 5): index of
                the first neutrino energy density column in self.endens
            useT (default False): if the run was performed considering
                a low reheating scenario,
                plot Neff as a function of t instead of x
        """
        if hasattr(self, "Neffdat") and not np.isnan(self.Neffdat[0, 0]):
            data = self.Neffdat
        else:
            try:
                rhogammas = self.endens[:, (endensx, endensgamma)][:, :]
            except (AttributeError, TypeError):
                print(traceback.format_exc())
                try:
                    self.zdat[:, (0, self.zCol)][:]
                except (AttributeError, TypeError):
                    print(traceback.format_exc())
                    return
                rhogammas = np.array(
                    [[x, PISQD15 * z ** 4] for x, z in self.zdat[:, (0, self.zCol)]]
                )
            try:
                rns = np.sum(
                    self.endens[:, endensnu0 : endensnu0 + self.nnu][:, :], axis=1
                )
                xs = self.endens[:, endensx][:]
                rhonus = np.array(list(zip(xs, rns)))
            except (AttributeError, TypeError):
                print(traceback.format_exc())
                try:
                    [self.rho[inu, inu, 0][:] for inu in range(self.nnu)]
                except (AttributeError, TypeError):
                    print(traceback.format_exc())
                    return
                rhonus = np.array(
                    [
                        [
                            x,
                            np.sum(
                                [
                                    self.integrateRho_yn(inu, 3, ix=ix)
                                    for inu in range(self.nnu)
                                ]
                            ),
                        ]
                        for ix, [x, z] in enumerate(self.zdat[:, (0, self.zCol)])
                    ]
                )
            frhonu = interp1d(*stripRepeated(rhonus, 0, 1))
            frhoga = interp1d(*stripRepeated(rhogammas, 0, 1))
            if self.lowReheating:
                try:
                    tv = self.endens[:, endensx + 1]
                except (AttributeError, TypeError):
                    print(traceback.format_exc())
                    try:
                        tv = self.zdat[:, 1]
                    except (AttributeError, TypeError):
                        print(traceback.format_exc())
                        return
                data = np.asarray(
                    [
                        [
                            x,
                            t,
                            8.0 / 7.0 * frhonu(x) / frhoga(x),
                            8.0
                            / 7.0
                            * frhonu(x)
                            / frhoga(x)
                            * (11.0 / 4.0) ** (4.0 / 3.0),
                        ]
                        for x, t in zip(rhogammas[:, 0], tv)
                    ]
                )
            else:
                data = np.asarray(
                    [
                        [
                            x,
                            8.0 / 7.0 * frhonu(x) / frhoga(x),
                            8.0
                            / 7.0
                            * frhonu(x)
                            / frhoga(x)
                            * (11.0 / 4.0) ** (4.0 / 3.0),
                        ]
                        for x in rhogammas[:, 0]
                    ]
                )
            print("Saving Neff to %s" % os.path.join(self.folder, "Neff.dat"))
            np.savetxt(os.path.join(self.folder, "Neff.dat"), data, fmt="%.7e")
        plt.plot(
            *stripRepeated(
                data,
                1 if self.lowReheating and useT else 0,
                2 if self.lowReheating else 1,
            ),
            ls=ls,
            c=lc,
            label=self.label if lab is None else lab
        )
        plt.xscale("log")
        plt.xlabel("$x$")
        if axes:
            ax = plt.gca()
            ax.set_ylim(nefflims)
            ax1 = ax.twinx()
            ax.set_ylabel(r"$N_{\rm eff}^{\rm in}$")
            ax1.set_ylabel(r"$N_{\rm eff}^{\rm now}$")
            ax1.set_ylim(np.asarray(nefflims) * (11.0 / 4) ** (4.0 / 3))

    def plotEnergyDensity(
        self,
        gamma_e=True,
        gec="#00ccff",
        ges=":",
        gamma_e_mu=True,
        gemc="#6666ff",
        gems="--",
        labels=[
            r"$\gamma$",
            "$e$",
            r"$\mu$",
            r"$\nu_e$",
            r"$\nu_\mu$",
            r"$\nu_\tau$",
            r"$\nu_s$",
        ],
        colors=["r", "b", "g", "#ff9933", "#ff9933", "#ff9933", "#ff00ff"],
        styles=["-", "-", "-", ":", "-.", "--", "-"],
        skip=[False, False, False, False, False, False, False],
        lw=1,
        allstyles=None,
        alllabels=None,
    ):
        """Plot the evolution of the single components
        and of the total energy density

        Parameters:
            gamma_e (default True): it True, plot the sum
                of the energies of photons and electrons
            gec (default "#00ccff"): color for the line of the sum
                of the energies of photons and electrons
            ges (default ":"): style for the line of the sum
                of the energies of photons and electrons
            gamma_e_mu (default True):it True, plot the sum
                of the energies of photons, electrons and muons
            gemc (default "#6666ff"): color for the line of the sum
                of the energies of photons, electrons and muons
            gems (default "--"): style for the line of the sum
                of the energies of photons, electrons and muons
            labels (default [
                    r"$\gamma$",
                    "$e$",
                    r"$\mu$",
                    r"$\nu_e$",
                    r"$\nu_\mu$",
                    r"$\nu_\tau$",
                    r"$\nu_s$",
                ]):
                the list of labels for all the lines
            colors (default
                ["r", "b", "g", "#ff9933", "#ff9933", "#ff9933", "#ff00ff"]):
                a list of colors for each line
            styles (default ["-", "-", "-", ":", "-.", "--", "-"]):
                a list of styles for each line
            skip (default [False]*7): True or False for each line
                to skip it and do not plot it
            lw (default 1): line width for the lines
            allstyles (default None): if it evaluates to True,
                a common line style for the all the lines
            alllabels (default None): if it evaluates to True,
                a common label for all the lines
        """
        try:
            self.endens[:, 0]
        except (AttributeError, TypeError):
            print(traceback.format_exc())
            return
        plt.plot(
            self.endens[:, 0],
            np.asarray([np.sum(cl[self.zCol + 1 :]) for cl in self.endens]),
            label="total" if alllabels is None else alllabels,
            c="k",
            ls="-" if not allstyles else allstyles,
            lw=lw,
        )
        if (
            self.lowReheating
            and self.endens.shape[1] == 10
            and labels[-1] == r"$\nu_s$"
        ):
            labels[-1] = r"$\phi$"
        for ix, lab in enumerate(labels):
            if skip[ix]:
                continue
            try:
                plt.plot(
                    self.endens[:, 0],
                    self.endens[:, self.zCol + 1 + ix],
                    label=lab if alllabels is None else alllabels,
                    c=colors[ix],
                    ls=styles[ix] if not allstyles else allstyles,
                    lw=lw,
                )
            except IndexError:
                pass
        if gamma_e:
            plt.plot(
                self.endens[:, 0],
                self.endens[:, self.zCol + 1] + self.endens[:, self.zCol + 2],
                label=r"$\gamma+e$" if alllabels is None else alllabels,
                c=gec,
                ls=ges if not allstyles else allstyles,
                lw=lw,
            )
        if gamma_e_mu:
            plt.plot(
                self.endens[:, 0],
                self.endens[:, self.zCol + 1]
                + self.endens[:, self.zCol + 2]
                + self.endens[:, self.zCol + 3],
                label=r"$\gamma+e+\mu$" if alllabels is None else alllabels,
                c=gemc,
                ls=gems if not allstyles else allstyles,
                lw=lw,
            )

    def plotEntropy(
        self,
        gamma_e=True,
        gec="#00ccff",
        ges=":",
        gamma_e_mu=True,
        gemc="#6666ff",
        gems="--",
        labels=[
            r"$\gamma$",
            "$e$",
            r"$\mu$",
            r"$\nu_e$",
            r"$\nu_\mu$",
            r"$\nu_\tau$",
            r"$\nu_s$",
        ],
        colors=["r", "b", "g", "#ff9933", "#ff9933", "#ff9933", "#ff00ff"],
        styles=["-", "-", "-", ":", "-.", "--", "-"],
        skip=[False, False, False, False, False, False, False],
        lw=1,
        allstyles=None,
        alllabels=None,
    ):
        """Plot the evolution of the single components
        and of the total entropy density

        Parameters:
            gamma_e (default True): it True, plot the sum
                of the entropies of photons and electrons
            gec (default "#00ccff"): color for the line of the sum
                of the entropies of photons and electrons
            ges (default ":"): style for the line of the sum
                of the entropies of photons and electrons
            gamma_e_mu (default True):it True, plot the sum
                of the entropies of photons, electrons and muons
            gemc (default "#6666ff"): color for the line of the sum
                of the entropies of photons, electrons and muons
            gems (default "--"): style for the line of the sum
                of the entropies of photons, electrons and muons
            labels (default [
                    r"$\gamma$",
                    "$e$",
                    r"$\mu$",
                    r"$\nu_e$",
                    r"$\nu_\mu$",
                    r"$\nu_\tau$",
                    r"$\nu_s$",
                ]):
                the list of labels for all the lines
            colors (default
                ["r", "b", "g", "#ff9933", "#ff9933", "#ff9933", "#ff00ff"]):
                a list of colors for each line
            styles (default ["-", "-", "-", ":", "-.", "--", "-"]):
                a list of styles for each line
            skip (default [False]*7): True or False for each line
                to skip it and do not plot it
            lw (default 1): line width for the lines
            allstyles (default None): if it evaluates to True,
                a common line style for the all the lines
            alllabels (default None): if it evaluates to True,
                a common label for all the lines
        """
        try:
            self.entropy[:, 0]
        except (AttributeError, TypeError):
            print(traceback.format_exc())
            return
        plt.plot(
            self.entropy[:, 0],
            np.asarray([np.sum(cl[self.zCol + 1 :]) for cl in self.entropy]),
            label="total" if alllabels is None else alllabels,
            c="k",
            ls="-" if not allstyles else allstyles,
            lw=lw,
        )
        for ix, lab in enumerate(labels):
            if skip[ix]:
                continue
            try:
                plt.plot(
                    self.entropy[:, 0],
                    self.entropy[:, self.zCol + 1 + ix],
                    label=lab if alllabels is None else alllabels,
                    c=colors[ix],
                    ls=styles[ix] if not allstyles else allstyles,
                    lw=lw,
                )
            except IndexError:
                pass
        if gamma_e:
            plt.plot(
                self.entropy[:, 0],
                self.entropy[:, self.zCol + 1] + self.entropy[:, self.zCol + 2],
                label=r"$\gamma+e$" if alllabels is None else alllabels,
                c=gec,
                ls=ges if not allstyles else allstyles,
                lw=lw,
            )
        if gamma_e_mu:
            plt.plot(
                self.entropy[:, 0],
                self.entropy[:, self.zCol + 1]
                + self.entropy[:, self.zCol + 2]
                + self.entropy[:, self.zCol + 3],
                label=r"$\gamma+e+\mu$" if alllabels is None else alllabels,
                c=gemc,
                ls=gems if not allstyles else allstyles,
                lw=lw,
            )

    def plotNumberDensity(
        self,
        labels=[
            r"$\gamma$",
            "$e$",
            r"$\mu$",
            r"$\nu_e$",
            r"$\nu_\mu$",
            r"$\nu_\tau$",
            r"$\nu_s$",
        ],
        colors=["r", "b", "g", "#ff9933", "#ff9933", "#ff9933", "#ff00ff"],
        styles=["-", "-", "-", ":", "-.", "--", "-"],
        skip=[False, False, False, False, False, False, False],
        lw=1,
        allstyles=None,
        alllabels=None,
    ):
        """Plot the evolution of the number density of each species

        Parameters:
            labels (default [
                    r"$\gamma$",
                    "$e$",
                    r"$\mu$",
                    r"$\nu_e$",
                    r"$\nu_\mu$",
                    r"$\nu_\tau$",
                    r"$\nu_s$",
                ]):
                the list of labels for all the lines
            colors (default
                ["r", "b", "g", "#ff9933", "#ff9933", "#ff9933", "#ff00ff"]):
                a list of colors for each line
            styles (default ["-", "-", "-", ":", "-.", "--", "-"]):
                a list of styles for each line
            skip (default [False]*7): True or False for each line
                to skip it and do not plot it
            lw (default 1): line width for the lines
            allstyles (default None): if it evaluates to True,
                a common line style for the all the lines
            alllabels (default None): if it evaluates to True,
                a common label for all the lines
        """
        try:
            self.number[:, 0]
        except (AttributeError, TypeError):
            print(traceback.format_exc())
            return
        for ix, lab in enumerate(labels):
            if skip[ix]:
                continue
            try:
                plt.plot(
                    self.number[:, 0],
                    self.number[:, self.zCol + 1 + ix],
                    label=lab if alllabels is None else alllabels,
                    c=colors[ix],
                    ls=styles[ix] if not allstyles else allstyles,
                    lw=lw,
                )
            except IndexError:
                pass

    def plotPArthENoPE(self, x="x", y="drhonu_dx", filter99=True, **kwargs):
        """Plot quantities that are employed in the PArthENoPE code for BBN

        Parameters:
            x (default "x"): the name of the quantity to plot in the x axis
            y (default "drhonu_dx"): the name of the quantity to plot in the y axis
            filter99 (default True): apply the filter to consider only
                points for which rho_rad > 0.99 rho_tot
            **kwargs: parameters to be passed to plt.plot
        """
        if not self.hasBBN:
            return
        allowed_x = ["x", "z", "Tgamma"]
        if self.lowReheating:
            allowed_x.append("t")
        allowed_y = ["rhonu", "drhonu_dx", "drhonu_dx_savgol", "N_func", "N_savgol"]
        for f, a in [[x, allowed_x], [y, allowed_y]]:
            if f not in a:
                raise AttributeError(
                    "Cannot plot the requested quantity '%s'. " % f
                    + "Allowed ones are: %s" % a
                )
        xv = getattr(self, x)[self.filter99] if filter99 else getattr(self, x)
        yv = getattr(self, y)[self.filter99] if filter99 else getattr(self, y)
        plt.plot(xv, yv, **kwargs)

    def doAllPlots(self, yref=5.0, color="k"):
        """Produce a series of plots for the given simulation

        Parameters:
            yref (default 5.0): the reference value of y
                to use in the plots
            color (default 'k'): color to use in the z and w plot
        """
        plt.close()
        self.plotZ(lc=color, lab="z")
        self.plotW(lc=color, ls=":", lab="w")
        finalizePlot("%s/z.pdf" % self.folder, xlab="$x$", ylab=r"$z$", xscale="log")

        for i in range(self.nnu):
            self.plotRhoDiagY(i, yref, styles[i], lc=colors[i])
        finalizePlot(
            "%s/rho_diag.pdf" % self.folder,
            xlab="$x$",
            ylab=r"$\rho_{\alpha\alpha}$",
            xscale="log",
            yscale="log",
        )

        for i in range(self.nnu):
            self.plotRhoDiagY(i, yref, styles[i], lc=colors[i], mass=True)
        finalizePlot(
            "%s/rho_mass_diag.pdf" % self.folder,
            xlab="$x$",
            ylab=r"$\rho_{\alpha\alpha}$",
            xscale="log",
            yscale="log",
        )

        for i in range(self.nnu):
            self.plotRhoFin(i, ls=styles[i], lc=colors[i], y2=True)
        finalizePlot("%s/rhofin_diag.pdf" % self.folder, xscale="linear", yscale="log")

        for i in range(self.nnu):
            self.plotRhoFin(i, ls=styles[i], lc=colors[i], y2=True, mass=True)
        finalizePlot(
            "%s/rhofin_mass_diag.pdf" % self.folder, xscale="linear", yscale="log"
        )

        if self.full:
            for i in range(self.nnu):
                for j in range(i + 1, self.nnu):
                    self.plotRhoOffDiagY(i, j, yref, lc=colors[2 * i + j - 1])
            finalizePlot("%s/rho_offdiag.pdf" % self.folder)
            for i in range(self.nnu):
                for j in range(i + 1, self.nnu):
                    self.plotdRhoOffDiagY(i, j, yref, lc=colors[2 * i + j - 1])
            finalizePlot("%s/drho_offdiag.pdf" % self.folder)

        self.plotNeff(lc=color, axes=False)
        finalizePlot("%s/Neff.pdf" % self.folder, legend=False, Neff_axes=True)

    def integrateRho_yn(self, inu, n, ix=-1, show=False, mass=False):
        """Compute the integral
        Int_0^Inf dy y^n f(y)/Pi^2
        for the requested eigenstate at the given x

        Parameters:
            inu: the neutrino eigenstate index
            n: the power to use in the exponential
            ix (default -1): the index of the requested x
                from the FortEPiaNO output
            show (default False): if True, print the result
                before returning it
            mass (default False): if True, use the density matrix
                in the mass basis instead of in the flavor one

        Output:
            the result of the integral
        """
        if mass:
            rho = self.rhoM
        else:
            rho = self.rho
        fy = interp1d(self.yv, rho[inu, inu, 0][ix, 1:] * (np.exp(self.yv) + 1))
        res = quad(lambda y: y ** n * fy(y) / (np.exp(y) + 1), self.yv[0], self.yv[-1])
        if show:
            print(res)
        return res[0] / np.pi ** 2
