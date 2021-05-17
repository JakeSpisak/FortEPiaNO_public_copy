"""Functions and classes that read the output and help to do plots"""

import os
import re
import sys
import traceback
import argparse

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
    import pandas as pd
except ImportError:
    print("Cannot import pandas...may raise errors later")
    pd = None
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
ELECTRONMASS_MEV = 0.51099895
muonLabel = r"$\mu$"


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
    inTicks=True,
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
    if not Neff_axes and inTicks:
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
        label="",
        deltas=False,
        full=True,
        plots=False,
        verbose=True,
    ):
        """Read the entire output of FortEPiaNO from a specific folder.
        It will ignore non-existing files and store all the available
        information for further processing (plots, ...)

        Parameters:
            folder: the name of the folder to consider
            nnu (default 3): number of neutrinos to consider
                (if more rows/columns of the density matrix exist,
                they will be ignored)
            label (default ""): a label to assign
                to the current FortEPiaNO point in the plots
            deltas (default False): if True, print the relative variation
                of energy and number density for each neutrino
            full (default True): if True, read also all the off-diagonal
                density matrix elements, otherwise ignore them
                (to save time if not needed in the plots, for example)
            plots (default False): if True, produce a series of plots
                after having read all the files
            verbose (default True): if True, print more error messages
                (e.g. when the folder is not found or w is not saved)
        """
        self.folder = folder
        self.full = full
        self.label = label
        self.verbose = verbose
        self.nnu = nnu
        self.ens_header = []
        if not os.path.exists(folder):
            if verbose:
                print("non-existing folder: %s" % folder)
            return
        self.readIni()
        self.readFD()
        self.readWZ()
        self.readNeff()
        self.readHeaders()
        self.readEENDensities(deltas)
        self.readResume()
        self.readNuDensMatrix(full)
        self.prepareRhoFinal()
        self.printTableLine()
        if plots:
            self.doAllPlots()

    def loadtxt(self, fname, *args, **kwargs):
        """Read a text file with numpy.loadtxt or pandas.read_csv"""
        try:
            return np.array(
                pd.read_csv(fname, *args, header=None, delimiter="\s+", **kwargs).values
            )
        except AttributeError:
            return np.loadtxt(fname, *args, **kwargs)

    def checkZdat(self):
        """Check if zdat has been read from the file or not.
        If it contains np.nan, assume self.zdat was not read from file
        """
        return (
            hasattr(self, "zdat")
            and isinstance(self.zdat, np.ndarray)
            and not np.isnan(self.zdat).all()
        )

    @property
    def x(self):
        """return x"""
        return self.zdat[:, 0] if self.checkZdat() else np.nan

    @property
    def z(self):
        """return z"""
        return self.zdat[:, 1] if self.checkZdat() else np.nan

    @property
    def zCol(self):
        """return the column index for z in the energy density file"""
        try:
            return self.ens_header.index("z")
        except ValueError:
            if self.verbose:
                print("z not in header file!")
            return 1

    @property
    def phCol(self):
        """return the column index for photons in the energy density file"""
        try:
            return self.ens_header.index("photon")
        except ValueError:
            if self.verbose:
                print("photon not in header file!")
            return self.zCol + 1

    @property
    def eCol(self):
        """return the column index for electrons in the energy density file"""
        try:
            return self.ens_header.index("electron")
        except ValueError:
            if self.verbose:
                print("electron not in header file!")
            return self.zCol + 2

    @property
    def muCol(self):
        """return the column index for muons in the energy density file"""
        try:
            return self.ens_header.index("muon")
        except ValueError:
            if self.verbose:
                print("muon not in header file!")
            return self.eCol + 1

    @property
    def nu1Col(self):
        """return the column index for nu1 in the energy density file"""
        try:
            return self.ens_header.index("nu1")
        except ValueError:
            if self.verbose:
                print("nu1 not in header file!")
            return self.zCol + 3

    @property
    def hasMuon(self):
        """return the column index for muons in the energy density file"""
        return "muon" in self.ens_header

    @property
    def Tgamma(self):
        """return T_gamma"""
        return self.z * ELECTRONMASS_MEV / self.x

    @property
    def t(self):
        """return t, if available"""
        return np.nan

    @property
    def rhonu(self):
        """return rhonu (total neutrino energy density)"""
        return np.nan

    @property
    def drhonu_dx(self):
        """return drhonu/dx"""
        return np.nan

    @property
    def drhonu_dx_savgol(self):
        """return a filtered version of drhonu/dx"""
        return np.nan

    @property
    def N_func(self):
        return np.nan

    @property
    def N_savgol(self):
        return np.nan

    def readEENDensities(self, deltas=False):
        """Read energy, entropy and number density,
        if the files are available

        Parameters:
            deltas (default False): if True, print the relative variation
                of energy and number density for each neutrino
        """
        refIx = 0
        try:
            self.endens = self.loadtxt("%s/energyDensity.dat" % self.folder)
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
                            self.endens[-1, self.nu1Col + i]
                            - self.endens[refIx, self.nu1Col + i]
                        )
                        / self.endens[refIx, self.nu1Col + i]
                        * 100
                        for i in range(self.nnu)
                    ]
                    self.tot_delta_ed = (
                        (
                            np.sum(
                                self.endens[-1, self.nu1Col : self.nu1Col + self.nnu]
                            )
                            - np.sum(
                                self.endens[refIx, self.nu1Col : self.nu1Col + self.nnu]
                            )
                        )
                        / np.sum(
                            self.endens[refIx, self.nu1Col : self.nu1Col + self.nnu]
                        )
                        * 100
                    )
                    if self.verbose:
                        print(
                            "delta energy density:\t"
                            + "\t".join(
                                [
                                    "nu%d: %f%%"
                                    % (
                                        i + 1,
                                        self.delta_ed[i],
                                    )
                                    for i in range(self.nnu)
                                ]
                            )
                        )
        try:
            self.entropy = self.loadtxt("%s/entropy.dat" % self.folder)
        except (IOError, OSError):
            self.entropy = np.nan
        else:
            try:
                self.entropy[:, :]
            except IndexError:
                pass
            else:
                if deltas:
                    ds = np.asarray([np.sum(cl[self.phCol :]) for cl in self.entropy])
                    self.tot_delta_sd = (ds[-1] - ds[refIx]) / ds[refIx] * 100
                    if self.verbose:
                        print("delta entropy density:\t%f%%" % self.tot_delta_sd)
        try:
            self.number = self.loadtxt("%s/numberDensity.dat" % self.folder)
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
                            self.number[-1, self.nu1Col + i]
                            - self.number[refIx, self.nu1Col + i]
                        )
                        / self.number[refIx, self.nu1Col + i]
                        * 100
                        for i in range(self.nnu)
                    ]
                    self.tot_delta_nd = (
                        (
                            np.sum(
                                self.number[-1, self.nu1Col : self.nu1Col + self.nnu]
                            )
                            - np.sum(
                                self.number[refIx, self.nu1Col : self.nu1Col + self.nnu]
                            )
                        )
                        / np.sum(
                            self.number[refIx, self.nu1Col : self.nu1Col + self.nnu]
                        )
                        * 100
                    )
                    if self.verbose:
                        print(
                            "delta number density:\t"
                            + "\t".join(
                                [
                                    "nu%d: %f%%" % (i + 1, self.delta_nd[i])
                                    for i in range(self.nnu)
                                ]
                            )
                        )

    def readFD(self):
        """Read and store the information from fd.dat"""
        try:
            fdy = self.loadtxt("%s/fd.dat" % self.folder)
        except (IOError, OSError):
            self.yv = np.nan
            self.fd = np.nan
        else:
            self.yv = fdy[:, 0]
            self.fd = fdy[:, 1]

    def readHeaders(self):
        """Read and store the information from the header files"""
        try:
            with open("%s/ens_header.dat" % self.folder) as _f:
                head = _f.read()
        except (IOError, OSError):
            print("Cannot read header. Assume the standard one")
            self.ens_header = "x z photon electron nu1 nu2 nu3".split()
        else:
            search = re.compile("nu \(1 to ([\d])\)")
            try:
                maxnu = int(search.findall(head)[0])
            except (AttributeError, IndexError):
                print("Cannot read number of neutrinos from header. Assume 3")
                maxnu = 3
            self.ens_header = head.replace(
                "nu (1 to %d)" % maxnu,
                " ".join(["nu%d" % (i + 1) for i in range(maxnu)]),
            ).split()

    def readIni(self):
        """Read and store the content of the ini.log file"""
        try:
            with open("%s/ini.log" % self.folder) as _ini:
                self.ini = _ini.read()
        except FileNotFoundError:
            self.ini = ""
        else:
            self.ini = self.ini.replace("\n", " ")
        # read flavorNumber
        search = re.compile("flavorNumber[\s]*=[\s]*([\d]+)")
        try:
            self.flavorNumber = int(search.findall(self.ini)[0])
        except (AttributeError, IndexError):
            print("Cannot read flavorNumber from ini.log. Assume 3")
            self.flavorNumber = 3
        if self.flavorNumber != self.nnu:
            print(
                "Warning: using nnu=%s, but FortEPiaNO used flavorNumber=%s"
                % (self.nnu, self.flavorNumber)
            )
        # read Ny
        search = re.compile("Ny[ ]*=[ ]*([\d]+)")
        try:
            self.Ny = int(search.findall(self.ini)[0])
        except (AttributeError, IndexError):
            print("Cannot read Ny from ini.log. Assume 3")
            self.Ny = 25

    def readIntermediate(self):
        """Read and store the intermediate quantities from the fortran code"""
        # read all the files and check that the number of lines are the same
        try:
            dat = self.loadtxt("%s/intermXF.dat" % self.folder)
        except (IOError, OSError):
            self.intermX = np.array([np.nan])
            self.intermN = np.array([np.nan])
        else:
            self.intermX = dat[:, 0]
            self.intermN = dat[:, 1]
        for f in [
            "intermY",
            "intermYdot",
            "intermHeff",
            "intermComm",
            "intermCTNue",
            "intermCTNunu",
        ]:
            try:
                setattr(self, f, self.loadtxt("%s/%s.dat" % (self.folder, f)))
            except (IOError, OSError):
                setattr(self, f, np.array([np.nan]))
            assert len(self.intermX) == len(getattr(self, f))

        try:
            self.nonRhoVars = self.intermY.shape[1] - self.intermHeff.shape[1]
        except IndexError:
            self.nonRhoVars = 0

        # separate quantities in Y and Ydot
        for t in ["", "dot"]:
            for f, ic in [
                ["intermZ", -1],
                ["intermW", -2],
            ]:
                try:
                    setattr(self, f + t, getattr(self, "intermY" + t)[:, ic])
                except IndexError:
                    setattr(self, f + t, np.array([np.nan]))
            try:
                setattr(
                    self,
                    "intermRho" + t,
                    getattr(self, "intermY" + t)[:, : -self.nonRhoVars],
                )
            except IndexError:
                setattr(self, "intermRho" + t, np.array([np.nan]))
        if len(self.intermX) == 1 and np.isnan(self.intermX[0]):
            return
        for f in [
            "intermRho",
            "intermRhodot",
            "intermHeff",
            "intermComm",
            "intermCTNue",
            "intermCTNunu",
        ]:
            setattr(
                self,
                f,
                np.array([self.reshapeVectorToMatrices(a) for a in getattr(self, f)]),
            )

    def readNeff(self):
        """Read the Neff.dat file"""
        try:
            self.Neffdat = self.loadtxt("%s/Neff.dat" % self.folder)
        except (IOError, OSError):
            self.Neffdat = np.asarray([[np.nan, np.nan, np.nan]])

    def readNuDensMatrix(self, full=True):
        """Read the evolution of the various components
        of the neutrino density matrix

        Parameters:
            full (default True): if True, read also all the off-diagonal
                density matrix elements, otherwise ignore them
        """
        self.rho = np.asarray(
            [[[None, None] for i in range(self.nnu)] for j in range(self.nnu)]
        )
        self.rhoM = np.asarray(
            [[[None, None] for i in range(self.nnu)] for j in range(self.nnu)]
        )
        for i in range(self.nnu):
            try:
                self.rho[i, i, 0] = self.loadtxt(
                    "%s/nuDens_diag%d.dat" % (self.folder, i + 1)
                )
            except (IOError, OSError):
                self.rho[i, i, 0] = np.nan
            if full:
                for j in range(i + 1, self.nnu):
                    try:
                        self.rho[i, j, 0] = self.loadtxt(
                            "%s/nuDens_nd_%d%d_re.dat" % (self.folder, i + 1, j + 1)
                        )
                    except (IOError, OSError):
                        self.rho[i, j, 0] = np.nan
                    try:
                        self.rho[i, j, 1] = self.loadtxt(
                            "%s/nuDens_nd_%d%d_im.dat" % (self.folder, i + 1, j + 1)
                        )
                    except (IOError, OSError):
                        self.rho[i, j, 1] = np.nan
            try:
                self.rhoM[i, i, 0] = self.loadtxt(
                    "%s/nuDens_mass%d.dat" % (self.folder, i + 1)
                )
            except (IOError, OSError):
                self.rhoM[i, i, 0] = np.nan

    def readResume(self):
        """Read the resume file and the information that it contains"""
        try:
            with open("%s/resume.dat" % self.folder) as _f:
                self.resume = _f.read()
        except FileNotFoundError:
            self.resume = ""
            self.hasResume = False
        else:
            self.resume = self.resume.replace("\n", " ")
            self.hasResume = True
        if self.hasResume:
            try:
                self.Neff = float(
                    re.search("Neff[ ]*=[ ]*([-\d.]*)", self.resume).group(1)
                )
            except (AttributeError, ValueError):
                self.Neff = np.nan
            try:
                self.wfin = float(
                    re.search("final w[ ]*=[ ]*([-\d.]*)", self.resume).group(1)
                )
            except (AttributeError, ValueError):
                if self.verbose:
                    print("cannot read w from resume.dat")
                self.wfin = np.nan
            try:
                self.zfin = float(
                    re.search("final z[ ]*=[ ]*([-\d.]*)", self.resume).group(1)
                )
            except (AttributeError, ValueError):
                self.zfin = np.nan
            self.deltaNeffi = np.array([np.nan for i in range(self.nnu)])
            if re.search("deltaNeff_([\d]+)[ ]*=[ ]*([-\d.]*)", self.resume):
                search = re.compile("deltaNeff_([\d]+)[ ]*=[ ]*([-\d.]*)")
                if search.findall(self.resume):
                    try:
                        self.deltaNeffi = np.array(
                            [float(g[1]) for g in search.findall(self.resume)]
                        )
                    except (AttributeError, IndexError):
                        if self.verbose:
                            print("cannot read deltaNeffi")
            else:
                # this part is kept for compatibility with previous versions
                search = re.compile("nuFactor([\d]+)[ ]*=[ ]*([E+\-\d.]*)")
                try:
                    factors = np.array([float(g[1]) for g in search.findall(self.ini)])
                except (AttributeError, IndexError):
                    factors = [np.nan] * 6
                deltarhofin = []
                search = re.compile("dRho_([\d]+)[ ]*=[ ]*([-\d.]*)")
                if search.findall(self.resume):
                    for i, g in enumerate(search.findall(self.resume)):
                        try:
                            deltarhofin.append(float(g[1]) + factors[i])
                        except (AttributeError, IndexError):
                            deltarhofin.append(np.nan)
                    deltarhofin = np.array(deltarhofin)
                    try:
                        self.deltaNeffi = self.Neff * deltarhofin / np.sum(deltarhofin)
                    except (AttributeError, ValueError):
                        if self.verbose:
                            print("cannot convert old dRho_i to new deltaNeffi")

    def readWZ(self):
        """Read the z.dat file"""
        try:
            self.zdat = self.loadtxt("%s/z.dat" % self.folder)
        except (IOError, OSError):
            self.zdat = np.asarray([[np.nan, np.nan, np.nan]])

    def reshapeVectorToMatrices(self, vec):
        """Take a Ny*flavorNumber or Ny*flavorNumber**2 vector
        and reshape it to a multidimensional numpy array,
        according to the order used by the fortran part of FortEPiaNO.
        The output array has shape (flavorNumber, flavorNumber, 2, Ny),
        where the first two dimensions identify the matrix element,
        the third dimension is for the real or imaginary part
        and the last dimension correspond to the momentum nodes.

        Parameter:
            vec: a list or 1D array,
                for the reshaping to work it must have length
                Ny*flavorNumber or Ny*flavorNumber**2

        Output:
            an array with shape (flavorNumber, flavorNumber, 2, Ny)
        """
        matrix = np.asarray(
            [
                [
                    [np.zeros(self.Ny), np.zeros(self.Ny)]
                    for i in range(self.flavorNumber)
                ]
                for j in range(self.flavorNumber)
            ]
        )
        if not isinstance(vec, np.ndarray):
            vec = np.array(vec)
        n = self.flavorNumber
        if len(vec) == self.Ny * self.flavorNumber:
            vec = vec.reshape(self.Ny, self.flavorNumber)
            for i in range(self.flavorNumber):
                matrix[i, i, 0] = vec[:, i]
        elif len(vec) == self.Ny * self.flavorNumber ** 2:
            vec = vec.reshape(self.Ny, self.flavorNumber ** 2)
            for i in range(self.flavorNumber):
                matrix[i, i, 0] = vec[:, i]
                for j in range(i + 1, self.flavorNumber):
                    vix = n + 2 * ((i * (n - 1) - int((i) * (i - 1) / 2)) + j - i - 1)
                    matrix[i, j, 0] = vec[:, vix]
                    matrix[i, j, 1] = vec[:, vix + 1]
                    matrix[j, i, 0] = matrix[i, j, 0]
                    matrix[j, i, 1] = -matrix[i, j, 1]
        else:
            raise ValueError("Invalid length of the array")
        return matrix

    def prepareRhoFinal(self):
        """Save the normalized diagonal entries of the final neutrino
        density matrix in mass basis, if not already existing
        """
        zid = (11.0 / 4.0) ** (1.0 / 3.0)

        def photonDensity(z):
            """photon energy density as a function of z"""
            return PISQD15 * z ** 4

        def Neffcontributions(delta, z, w, nnu, ymin=0.01, ymax=20):
            """neutrino energy density for each of the nnu eigenstates.
            Compute the integral (between ymin and ymax)
            of the final neutrino density matrix
            (diagonal entries) as a function of w,
            divide by the photon density (function of z)
            and adjust the normalization constants.
            The output
            """
            return [
                quad(
                    lambda y: y ** 3
                    / (np.exp(y / w) + 1)
                    * interp1d(delta[:, 0], delta[:, i], fill_value="extrapolate")(y),
                    ymin,
                    ymax,
                )[0]
                / np.pi ** 2
                * (zid) ** 4
                / photonDensity(z)
                / 0.875
                for i in range(1, nnu + 1)
            ]

        if not hasattr(self, "wfin") or np.isnan(self.wfin):
            if self.verbose:
                print("cannot read final w, I will not renormalize the final rho")
            return
        if not hasattr(self, "zfin") or np.isnan(self.zfin):
            if self.verbose:
                print("cannot read final z, I will not renormalize the final rho")
            return
        for fm in ["rho_final", "rho_final_mass"]:
            if os.path.exists("%s/%s.dat" % (self.folder, fm)) and not (
                os.path.exists("%s/%s_norm.dat" % (self.folder, fm))
                and os.path.exists("%s/%s_var.dat" % (self.folder, fm))
            ):
                data = self.loadtxt("%s/%s.dat" % (self.folder, fm))
                # Compute the variation of the neutrino density matrix
                # with respect to the FermiDirac at temperature w.
                # To obtain the contributions to Neff, integrate:
                # rhonu_i = 1/pi**2 \int_ymin^ymax dy y**3 var(y) / [exp(y/w)+1]
                # and normalize:
                # deltaNeff_i = rhonu_i * zid**4 / photonDensity(z) * 8/7
                # Here in this function, you can use:
                # deltaNeff_i = Neffcontributions(var, self.zfin, self.wfin, ...)
                var = np.column_stack(
                    (
                        data[:, 0],
                        np.array(
                            [
                                data[:, i] * (np.exp(data[:, 0] / self.wfin) + 1)
                                for i in range(1, data.shape[1])
                            ]
                        ).T,
                    )
                )
                # Compute the renormalized variation of the neutrino density matrix
                # with respect to the FermiDirac 1/[exp(y)+1].
                # To obtain the contributions to Neff, integrate:
                # rhonu_i = 1/pi**2 \int_ymin^ymax dy y**3 nor(y) / [exp(y)+1]
                # and normalize:
                # deltaNeff_i = rhonu_i * zid**4 / photonDensity(zid) * 8/7
                # Here in this function, you can use:
                # deltaNeff_i = Neffcontributions(nor, zid, 1., ...)
                nor = np.column_stack(
                    (
                        var[:, 0],
                        np.array(
                            [
                                interp1d(
                                    var[:, 0], var[:, i], fill_value="extrapolate"
                                )(var[:, 0] * self.wfin)
                                * (zid / self.zfin * self.wfin) ** 4
                                for i in range(1, data.shape[1])
                            ]
                        ).T,
                    )
                )
                try:
                    np.savetxt(
                        "%s/%s_var.dat" % (self.folder, fm),
                        var,
                        fmt="%15.7e",
                    )
                    np.savetxt(
                        "%s/%s_norm.dat" % (self.folder, fm),
                        nor,
                        fmt="%15.7e",
                    )
                except IOError:
                    print("Cannot write the converted neutrino density matrices!")

    def interpolateRhoIJ(self, i1, i2, y, ri=0, yexp=0, mass=False):
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
            yexp (default 0): if >0, multiply the output by y**yexp
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
            cy = fy(y) * (y ** yexp if yexp > 0 else 1.0)
            if cy != prevy:
                prevy = cy
                yv.append(prevy)
                xv.append(x)
        if x != xv[-1] or cy != yv[-1]:
            xv.append(x)
            yv.append(cy)
        return np.asarray(xv), np.asarray(yv)

    def interpolateRhoIJ_x(self, i1, i2, x, ri=0, yexp=0, mass=False):
        """Interpolate any entry of the density matrix at a given x,
        and return its value for all the y grid points

        Parameters:
            i1: row index of the density matrix to consider
            i2: column index of the density matrix to consider
            x: the x value at which to interpolate
            ri (default 0): it 0, use real part, if 1 the imaginary one
                (only for off-diagonal entries of the density matrix)
            yexp (default 0): if >0, multiply the output by y**yexp
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
                rho[i1, i2, ri][:, i + 1] * (y ** yexp if yexp > 0 else 1.0),
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
            [self.deltaNeffi[i] for i in range(self.nnu)]
        except (AttributeError, IndexError):
            print(traceback.format_exc())
            return
        if not self.verbose:
            return
        if self.hasResume:
            deltastr = ""
            for i in range(self.nnu):
                deltastr += "{:.5f} & ".format(self.deltaNeffi[i])
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

    def plotRhoDiag(self, inu, iy, ls, lc="k", lab=None, mass=False):
        """Plot one diagonal element of the density matrix
        at a given point (index iy) in the momentum grid,
        as a function of x

        Parameters:
            inu: diagonal index of the density matrix entry to consider
            iy: the index of the requested momentum in the momentum grid
            ls: the line style
            lc (default "k"): the line color
            lab (default None): if not None, the custom value
                for the line label
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
            label=r"%s $\alpha$=%d" % (self.label, inu + 1) if lab is None else lab,
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
        self,
        i1,
        i2=None,
        ri=0,
        ls="-",
        lc="k",
        yexp=0,
        lab=None,
        mass=False,
        divide_fd=False,
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
            yexp (default 0): if >0,
                multiply the diagonal elements times y**yexp
            lab (default None): if not None, the line label
            mass (default False): if True, use the density matrix
                in the mass basis
            divide_fd (default False): if True,
                divide by the Fermi Dirac (self.fd/self.yv**2)
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
        fyv = (
            (self.yv ** yexp * rho[i1, i2, ri][-1, 1:])
            if yexp > 0
            else rho[i1, i2, ri][-1, 1:]
        )
        plt.plot(
            self.yv,
            fyv / ((self.fd / self.yv ** 2) if divide_fd else 1.0),
            ls=ls,
            c=lc,
            label=label,
        )
        plt.xlabel("$y$")
        plt.ylabel(
            r"$%s\rho_{\alpha\beta}^{\rm fin}(y)$"
            % ("y^{%d}" % yexp if yexp > 0 else "")
        )

    def plotRhoX(
        self,
        i1,
        x,
        i2=None,
        ri=0,
        ls="-",
        lc="k",
        lab="",
        yexp=0,
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
            yexp (default 0): if >0,
                multiply the diagonal elements times y**yexp
            mass (default False): if True, use the density matrix
                in the mass basis
            divide_by (default 1.0): divide the rho values by the given
                float or array
            divide_fd (default False): if True,
                divide by the Fermi Dirac (self.fd/self.yv**2)
        """
        if i2 is None:
            i2 = i1
        if ri not in [0, 1]:
            ri = 0
        try:
            interp = self.interpolateRhoIJ_x(i1, i2, x, ri, yexp=yexp, mass=mass)
        except (AttributeError, TypeError):
            print(traceback.format_exc())
            return
        xv, yv = interp
        plt.plot(
            xv,
            yv / divide_by / ((self.fd / self.yv ** 2) if divide_fd else 1.0),
            ls=ls,
            c=lc,
            label=r"%s $\alpha\beta$=%d%d %s x=%f"
            % (self.label, i1 + 1, i2 + 1, "re" if ri == 0 else "im", x)
            if not lab
            else lab,
        )
        plt.xlabel("$y$")
        plt.ylabel(r"$%s\rho_{\alpha\beta}(y)$" % ("y^{%d}" % yexp if yexp > 0 else ""))

    def plotRhoDiagY(
        self, inu, y, ls, lc="k", lab=None, yexp=0, mass=False, divide_by=1.0
    ):
        """Plot one diagonal element of the density matrix at a given y
        as a function of x

        Parameters:
            inu: diagonal element of the density matrix to consider
            y: the y value at which the density matrix is interpolated
            ls: the line style
            lc (default "k"): the line color
            lab (default None): if not None, the line label
            yexp (default 0): if >0,
                multiply the diagonal elements times y**yexp
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
            np.asarray(yv) * (y ** yexp if yexp > 0 else 1.0) / divide_by,
            label=label,
            ls=ls,
            c=lc,
        )
        plt.xscale("log")
        plt.xlabel("$x$")
        plt.ylabel(r"$%s\rho_{\alpha\alpha}$" % ("y^{%d}" % yexp if yexp > 0 else ""))

    def plotdRhoDiagY(self, inu, y, ls, lc="k", lab=None, yexp=0, mass=False):
        """Plot the x derivative (np.gradient)
        of one diagonal element of the density matrix at a given y
        as a function of x

        Parameters:
            inu: diagonal element of the density matrix to consider
            y: the y value at which the density matrix is interpolated
            ls: the line style
            lc (default "k"): the line color
            lab (default None): if not None, the line label
            yexp (default 0): if >0,
                multiply the diagonal elements times y**yexp
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
            np.gradient(np.asarray(yv) * (y ** yexp if yexp > 0 else 1.0), x),
            label=label,
            ls=ls,
            c=lc,
        )
        plt.xscale("log")
        plt.xlabel("$x$")
        plt.ylabel(
            r"$d%s\rho_{\alpha\alpha}/dx$" % ("y^{%d}" % yexp if yexp > 0 else "")
        )

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
            data = np.asarray(
                [
                    [
                        x,
                        8.0 / 7.0 * frhonu(x) / frhoga(x),
                        8.0 / 7.0 * frhonu(x) / frhoga(x) * (11.0 / 4.0) ** (4.0 / 3.0),
                    ]
                    for x in rhogammas[:, 0]
                ]
            )
            print("Saving Neff to %s" % os.path.join(self.folder, "Neff.dat"))
            np.savetxt(os.path.join(self.folder, "Neff.dat"), data, fmt="%.7e")
        plt.plot(
            *stripRepeated(
                data,
                0,
                1,
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

    def plotNeffAtAllX(
        self,
        lc="k",
        ls="-",
        lab=None,
        endensx=0,
        endensgamma=2,
        endensnu0=5,
    ):
        """Plot the evolution of Neff as a function of x, computing Neff
        at each moment using rho_nu, rho_gamma, z, w.
        May not give perfectly correct results due to numerical
        instability in the calculation of z/w.

        Parameters:
            lc (default "k"): the line color
            ls (default "-"): the line style
            lab (default None): if not None, the line label
            endensx (default 0): index of
                the x column in self.endens
            endensgamma (default 2): index of
                the photon energy density column in self.endens
            endensnu0 (default 5): index of
                the first neutrino energy density column in self.endens
        """
        try:
            rhogammas = self.endens[:, (endensx, endensgamma)][:, :]
        except (AttributeError, TypeError, IndexError):
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
            rns = np.sum(self.endens[:, endensnu0 : endensnu0 + self.nnu][:, :], axis=1)
            xs = self.endens[:, endensx][:]
            rhonus = np.array(list(zip(xs, rns)))
        except (AttributeError, TypeError, IndexError):
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
        try:
            zs = self.zdat[:, (0, self.zCol)][:]
            ws = self.zdat[:, self.zCol + 1]
        except (AttributeError, TypeError, IndexError):
            print(traceback.format_exc())
            return
        frhonu = interp1d(*stripRepeated(rhonus, 0, 1))
        frhoga = interp1d(*stripRepeated(rhogammas, 0, 1))
        zs[:, 1] = zs[:, 1] / ws[:]
        fzow = interp1d(*stripRepeated(zs, 0, 1))
        data = np.asarray(
            [
                [
                    x,
                    8.0 / 7.0 * frhonu(x) / frhoga(x) * fzow(x) ** 4.0,
                ]
                for x in rhogammas[:, 0]
            ]
        )
        plt.plot(
            data[:, 0],
            data[:, 1],
            ls=ls,
            c=lc,
            label=self.label if lab is None else lab,
        )
        plt.xscale("log")
        plt.xlabel(r"$x$")
        plt.ylabel(r"$N_{\rm eff}$")

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
            muonLabel,
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
        if not self.hasMuon and muonLabel in labels:
            del labels[labels.index(muonLabel)]
        plt.plot(
            self.endens[:, 0],
            np.asarray([np.sum(cl[self.phCol :]) for cl in self.endens]),
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
                    self.endens[:, 0],
                    self.endens[:, self.phCol + ix],
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
                self.endens[:, self.phCol] + self.endens[:, self.eCol],
                label=r"$\gamma+e$" if alllabels is None else alllabels,
                c=gec,
                ls=ges if not allstyles else allstyles,
                lw=lw,
            )
        if gamma_e_mu and self.hasMuon:
            plt.plot(
                self.endens[:, 0],
                self.endens[:, self.phCol]
                + self.endens[:, self.eCol]
                + self.endens[:, self.muCol],
                label=r"$\gamma+e+\mu$" if alllabels is None else alllabels,
                c=gemc,
                ls=gems if not allstyles else allstyles,
                lw=lw,
            )

    def plotDeltaEntropy(
        self,
        lc="k",
        ls="-",
        lab=None,
    ):
        """Plot the variation in total entropy as a function of x

        Parameters:
            lc (default "k"): the line color
            ls (default "-"): the line style
            lab (default None): if not None, the line label
        """
        ds = np.asarray([np.sum(cl[self.phCol :]) for cl in self.entropy])
        plt.plot(
            self.entropy[:, 0],
            (ds / ds[0] - 1.0) * 100.0,
            c=lc,
            ls=ls,
            label=self.label if lab is None else lab,
        )
        plt.xscale("log")
        plt.xlabel("$x$")
        plt.ylabel(r"$\delta s_{\rm tot}$ [%]")

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
            muonLabel,
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
        if not self.hasMuon and muonLabel in labels:
            del labels[labels.index(muonLabel)]
        plt.plot(
            self.entropy[:, 0],
            np.asarray([np.sum(cl[self.phCol :]) for cl in self.entropy]),
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
                    self.entropy[:, self.phCol + ix],
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
                self.entropy[:, self.phCol] + self.entropy[:, self.eCol],
                label=r"$\gamma+e$" if alllabels is None else alllabels,
                c=gec,
                ls=ges if not allstyles else allstyles,
                lw=lw,
            )
        if gamma_e_mu and self.hasMuon:
            plt.plot(
                self.entropy[:, 0],
                self.entropy[:, self.phCol]
                + self.entropy[:, self.eCol]
                + self.entropy[:, self.muCol],
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
            muonLabel,
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
        if not self.hasMuon and muonLabel in labels:
            del labels[labels.index(muonLabel)]
        for ix, lab in enumerate(labels):
            if skip[ix]:
                continue
            try:
                plt.plot(
                    self.number[:, 0],
                    self.number[:, self.phCol + ix],
                    label=lab if alllabels is None else alllabels,
                    c=colors[ix],
                    ls=styles[ix] if not allstyles else allstyles,
                    lw=lw,
                )
            except IndexError:
                pass

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
            yscale="log",
        )
        for i in range(self.nnu):
            self.plotdRhoDiagY(i, yref, styles[i], lc=colors[i])
        finalizePlot(
            "%s/drho_diag.pdf" % self.folder,
        )

        for i in range(self.nnu):
            self.plotRhoDiagY(i, yref, styles[i], lc=colors[i], mass=True)
        finalizePlot(
            "%s/rho_mass_diag.pdf" % self.folder,
            yscale="log",
        )
        for i in range(self.nnu):
            self.plotdRhoDiagY(i, yref, styles[i], lc=colors[i], mass=True)
        finalizePlot(
            "%s/drho_mass_diag.pdf" % self.folder,
        )

        for i in range(self.nnu):
            self.plotRhoFin(i, ls=styles[i], lc=colors[i], yexp=2)
        finalizePlot("%s/rhofin_diag.pdf" % self.folder, xscale="linear", yscale="log")

        for i in range(self.nnu):
            self.plotRhoFin(i, ls=styles[i], lc=colors[i], yexp=2, mass=True)
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


def setParser():
    """Prepare the parser for reading the command line arguments

    Output:
        the parser
    """
    parser = argparse.ArgumentParser(prog="fortepianoOutput.py")
    parser.add_argument(
        "folder",
        help="the name of the output folder that contains"
        + " the results of the Fortran code",
    )
    parser.add_argument(
        "--nnu",
        type=int,
        default=3,
        help="number of neutrinos to consider"
        + " (if more rows/columns of the density matrix exist, "
        + "they will be ignored)",
    )
    parser.add_argument(
        "--label",
        default="",
        help="a label to identify the run in the plot legends",
    )
    parser.add_argument(
        "--deltas",
        action="store_true",
        help="if True, print the relative variation of "
        + " energy and number density for each neutrino",
    )
    parser.add_argument(
        "--full",
        action="store_false",
        help="if True, read also all the off-diagonal"
        + " density matrix elements, otherwise ignore them"
        + " (to save time if not needed in the plots, for example)",
    )
    parser.add_argument(
        "--plots",
        action="store_true",
        help="if True, produce a series of plots after having read all the files",
    )
    parser.add_argument(
        "--verbose",
        action="store_false",
        help="increase the number of messages printed by the code",
    )
    return parser


if __name__ == "__main__":
    parser = setParser()
    args = parser.parse_args(sys.argv[1:])
    FortEPiaNORun(
        args.folder,
        **{
            k: getattr(args, k)
            for k in ["nnu", "label", "deltas", "full", "plots", "verbose"]
        }
    )
