#!/usr/bin/python
import argparse
import os
import shutil
import sys
import numpy as np
import matplotlib

matplotlib.use("agg")
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

if sys.version_info[0] < 3:
    import unittest2 as unittest
    from mock import call, patch, MagicMock

    USE_AUTOSPEC_CLASS = False
else:
    import unittest
    from unittest.mock import call, patch, MagicMock

    USE_AUTOSPEC_CLASS = True


import fortepianoOutput as fpom
import prepareIni as pim

try:
    FileNotFoundError
except NameError:
    FileNotFoundError = IOError

allPlots = True


class Namespace:
    """Produce a namespace which has attributesto be used
    instead of dictionary keys,
    to mimic the output of parse_args.Namespace
    """

    def __init__(self, **kwargs):
        """Use the input dictionary to update the instance __dict__

        Parameters:
            any keyword arguments
        """
        self.__dict__.update(kwargs)


def setUpModule():
    """Generate a partially empty resume file for tests"""
    folder = "output/no/"
    if os.path.exists(folder):
        shutil.rmtree(folder)
    os.makedirs(folder)
    with open("%s/resume.dat" % folder, "w") as _f:
        _f.write("final w =  NaN\nfinal z =  1.5\nNeff    =  \n")


def tearDownModule():
    """Delete test folder"""
    for folder in [
        "output/no/",
        "output/nonexistent",
        "output/nonexistent1",
        "output/nonexistent2",
    ]:
        if os.path.exists(folder):
            shutil.rmtree(folder)


class FPTestCase(unittest.TestCase):
    """additional test functions"""

    def assertEqualArray(self, a, b):
        # type: (np.ndarray, np.ndarray) -> None
        """Assert that two np.ndarrays (a, b) are equal
        using np.all(a == b)

        Parameters:
            a, b: the two np.ndarrays to compare
        """
        self.assertTrue(np.allclose(a, b, equal_nan=True))


class TestFortepianoOutput(FPTestCase):
    """Testing some content of the FortepianoOutput module"""

    def test_attributes(self):
        """test attributes"""
        self.assertIsInstance(fpom.colors, list)
        self.assertIsInstance(fpom.styles, list)
        self.assertIsInstance(fpom.markers, list)
        self.assertIsInstance(fpom.PISQD15, float)

    def test_finalizePlot(self):
        """test finalizePlot"""
        plt.figure()
        ax = plt.gca()
        ax1 = ax.twiny()
        with patch("matplotlib.pyplot.title") as _tit, patch(
            "matplotlib.pyplot.gca", return_value=ax
        ) as _gca, patch("matplotlib.pyplot.Axes.tick_params") as _tpa, patch(
            "matplotlib.pyplot.legend"
        ) as _leg, patch(
            "matplotlib.pyplot.xlabel"
        ) as _xla, patch(
            "matplotlib.pyplot.ylabel"
        ) as _yla, patch(
            "matplotlib.pyplot.xscale"
        ) as _xsc, patch(
            "matplotlib.pyplot.yscale"
        ) as _ysc, patch(
            "matplotlib.pyplot.xlim"
        ) as _xli, patch(
            "matplotlib.pyplot.ylim"
        ) as _yli, patch(
            "matplotlib.pyplot.tight_layout"
        ) as _tig, patch(
            "matplotlib.pyplot.savefig", side_effect=[None, None, FileNotFoundError]
        ) as _sfi, patch(
            "matplotlib.pyplot.close"
        ) as _clo:
            fpom.finalizePlot("fname", legend=False)
            _tit.assert_called_once_with("")
            _gca.assert_called_once_with()
            _tpa.assert_called_once_with(
                "both",
                which="both",
                direction="in",
                left=True,
                right=True,
                top=True,
                bottom=True,
            )
            self.assertEqual(_leg.call_count, 0)
            self.assertEqual(_xla.call_count, 0)
            self.assertEqual(_yla.call_count, 0)
            self.assertEqual(_xsc.call_count, 0)
            self.assertEqual(_ysc.call_count, 0)
            self.assertEqual(_xli.call_count, 0)
            self.assertEqual(_yli.call_count, 0)
            _tig.assert_called_once_with(rect=(-0.025, -0.03, 1.015, 1.03))
            _sfi.assert_called_once_with("fname")
            _clo.assert_called_once_with()
            fpom.finalizePlot(
                "fname",
                title="tit",
                xlab="x",
                ylab="y",
                xscale="log",
                yscale="log",
                xlim=[0, 1],
                ylim=[10, 11],
                tightrect=(0, 1, 2, 3),
            )
            _tit.assert_any_call("tit")
            _leg.assert_called_once_with(loc="best", ncol=1)
            _xla.assert_called_once_with("x")
            _yla.assert_called_once_with("y")
            _xsc.assert_called_once_with("log")
            _ysc.assert_called_once_with("log")
            _xli.assert_called_once_with([0, 1])
            _yli.assert_called_once_with([10, 11])
            _tig.assert_any_call(rect=(0, 1, 2, 3))
            self.assertEqual(_tpa.call_count, 2)
            fpom.finalizePlot(
                "fname",
                legcol=3,
                lloc="upper right",
            )
            _leg.assert_any_call(loc="upper right", ncol=3)

        # x_T
        lims = [5, 10]
        with patch("matplotlib.pyplot.title") as _tit, patch(
            "matplotlib.pyplot.gca", return_value=ax
        ) as _gca, patch("matplotlib.pyplot.Axes.tick_params") as _tpa, patch(
            "matplotlib.pyplot.Axes.get_xlim", return_value=lims
        ) as _gxl, patch(
            "matplotlib.pyplot.Axes.set_xscale"
        ) as _sxs, patch(
            "matplotlib.pyplot.Axes.set_xlabel"
        ) as _sxb, patch(
            "matplotlib.pyplot.Axes.twiny", return_value=ax1
        ) as _twy, patch(
            "matplotlib.pyplot.Axes.set_xlim"
        ) as _sxi, patch(
            "matplotlib.pyplot.tight_layout"
        ) as _tig, patch(
            "matplotlib.pyplot.savefig", side_effect=[None, None, FileNotFoundError]
        ) as _sfi, patch(
            "matplotlib.pyplot.close"
        ) as _clo:
            fpom.finalizePlot("fname", x_T=True)
            _tpa.assert_any_call(
                "both",
                which="both",
                direction="in",
                left=True,
                right=True,
                top=True,
                bottom=True,
            )
            _gxl.assert_called_once_with()
            _twy.assert_called_once_with()
            _sxb.assert_has_calls([call("$x$"), call("$T$ [MeV]")])
            _sxs.assert_has_calls([call("log"), call("log")])
            _sxi.assert_called_once_with(
                [0.5109989461 / lims[0], 0.5109989461 / lims[1]]
            )

        # Neff_axes
        with patch("matplotlib.pyplot.title") as _tit, patch(
            "matplotlib.pyplot.gca", return_value=ax
        ) as _gca, patch("matplotlib.pyplot.Axes.tick_params") as _tpa, patch(
            "matplotlib.pyplot.Axes.set_ylabel"
        ) as _syb, patch(
            "matplotlib.pyplot.Axes.get_ylim", return_value=lims
        ) as _gyi, patch(
            "matplotlib.pyplot.Axes.twinx", return_value=ax1
        ) as _twx, patch(
            "matplotlib.pyplot.Axes.set_ylim"
        ) as _syi, patch(
            "matplotlib.pyplot.tight_layout"
        ) as _tig, patch(
            "matplotlib.pyplot.savefig", side_effect=[None, None, FileNotFoundError]
        ) as _sfi, patch(
            "matplotlib.pyplot.close"
        ) as _clo:
            fpom.finalizePlot("fname", Neff_axes=True)
            _gyi.assert_called_once_with()
            _twx.assert_called_once_with()
            self.assertEqual(_tpa.call_count, 2)
            _tpa.assert_has_calls(
                [
                    call(
                        "both",
                        which="both",
                        direction="in",
                        left=True,
                        right=False,
                        labelleft=True,
                        labelright=False,
                    ),
                    call(
                        "both",
                        which="both",
                        direction="in",
                        left=False,
                        right=True,
                        labelleft=False,
                        labelright=True,
                    ),
                ]
            )
            _syb.assert_has_calls(
                [
                    call(
                        r"$N_{\rm eff}^{\rm in}=\frac{8}{7}\frac{\rho_\nu}{\rho_\gamma}$"
                    ),
                    call(
                        r"$N_{\rm eff}^{\rm now}="
                        + r"\frac{8}{7}\left(\frac{11}{4}\right)^{4/3}\;"
                        + r"\frac{\rho_\nu}{\rho_\gamma}$"
                    ),
                ]
            )
            self.assertEqualArray(
                _syi.call_args[0][0], np.asarray(lims) * (11.0 / 4) ** (4.0 / 3)
            )

    def test_stripRepeated(self):
        """test stripRepeated"""
        data = [[0, 1, 2, 3], [10, 11, 12, 13], [20, 11, 22, 23], [30, 31, 22, 33]]
        x, y = fpom.stripRepeated(data, 0, 1)
        self.assertIsInstance(x, np.ndarray)
        self.assertIsInstance(y, np.ndarray)
        self.assertEqualArray(x, [0, 10, 30])
        self.assertEqualArray(y, [1, 11, 31])
        x, y = fpom.stripRepeated(data, 3, 1)
        self.assertEqualArray(x, [3, 13, 33])
        self.assertEqualArray(y, [1, 11, 31])
        x, y = fpom.stripRepeated(data, 3, 2)
        self.assertEqualArray(x, [3, 13, 23, 33])
        self.assertEqualArray(y, [2, 12, 22, 22])
        x, y = fpom.stripRepeated(data, 0, 4)
        self.assertIsInstance(x, np.ndarray)
        self.assertIsInstance(y, np.ndarray)
        self.assertEqualArray(x, [np.nan])
        self.assertEqualArray(y, [np.nan])

    def test_setParser(self):
        """test the fortepianoOutput.setParser function"""
        parser = argparse.ArgumentParser()
        parser.add_argument = MagicMock()
        with patch("argparse.ArgumentParser", return_value=parser) as _ap, patch(
            "argparse.ArgumentParser.add_argument"
        ) as _aa:
            self.assertEqual(fpom.setParser(), parser)
            _ap.assert_called_once_with(prog="fortepianoOutput.py")
        parser.add_argument.assert_has_calls(
            [
                call(
                    "folder",
                    help="the name of the output folder that contains"
                    + " the results of the Fortran code",
                ),
                call(
                    "--nnu",
                    type=int,
                    default=3,
                    help="number of neutrinos to consider"
                    + " (if more rows/columns of the density matrix exist, "
                    + "they will be ignored)",
                ),
                call(
                    "--label",
                    default="",
                    help="a label to identify the run in the plot legends",
                ),
                call(
                    "--deltas",
                    action="store_true",
                    help="if True, print the relative variation of "
                    + " energy and number density for each neutrino",
                ),
                call(
                    "--full",
                    action="store_false",
                    help="if True, read also all the off-diagonal"
                    + " density matrix elements, otherwise ignore them"
                    + " (to save time if not needed in the plots, for example)",
                ),
                call(
                    "--plots",
                    action="store_true",
                    help="if True, produce a series of plots"
                    + " after having read all the files",
                ),
                call(
                    "--verbose",
                    action="store_false",
                    help="increase the number of messages printed by the code",
                ),
            ],
            any_order=True,
        )

    def test_parsing(self):
        """test that there are no errors in the parsing process"""
        parser = fpom.setParser()
        with self.assertRaises(SystemExit):
            args = parser.parse_args(["abc", "def"])
        args = parser.parse_args(["outdir"])
        self.assertEqual(args.folder, "outdir")
        self.assertEqual(args.nnu, 3)
        self.assertEqual(args.label, "")
        self.assertEqual(args.deltas, False)
        self.assertEqual(args.full, True)
        self.assertEqual(args.plots, False)
        self.assertEqual(args.verbose, True)
        for arg, val, lin in [
            ["nnu", 4, "--nnu=4"],
            ["label", "'abc d'", "--label='abc d'"],
            ["deltas", True, "--deltas"],
            ["full", False, "--full"],
            ["plots", True, "--plots"],
            ["verbose", False, "--verbose"],
        ]:
            args = parser.parse_args(["outdir", lin])
            self.assertEqual(args.folder, "outdir")
            self.assertEqual(args.nnu, val if "nnu" == arg else 3)
            self.assertEqual(args.label, val if "label" == arg else "")
            self.assertEqual(args.deltas, val if "deltas" == arg else False)
            self.assertEqual(args.full, val if "full" == arg else True)
            self.assertEqual(args.plots, val if "plots" == arg else False)
            self.assertEqual(args.verbose, val if "verbose" == arg else True)
        for lin in [
            "--nnu=4.25",
            "--nnu=aaa",
        ]:
            with self.assertRaises(SystemExit):
                args = parser.parse_args(["outdir", lin])


class TestFortEPiaNORun(FPTestCase):
    """Testing the FortEPiaNORun class"""

    @classmethod
    def setUpClass(self):
        """Set maxDiff to None and load the output of explanatory.ini"""
        self.maxDiff = None
        self.emptyRun = fpom.FortEPiaNORun("/nonexistent")
        self.explanatory = fpom.FortEPiaNORun("output/", label="label")
        self.failedRun = fpom.FortEPiaNORun("output/no/")

    def runAllPlots(self, run):
        """Call all the plotting function from FortEPiaNORun"""
        if not allPlots:
            return
        run.plotFD()
        run.doAllPlots()
        run.plotZoverW()
        run.plotDeltaZ(run)
        run.plotRhoDiag(0, 4, "-")
        run.plotdRhoDiag(0, 4, "-")
        run.plotRhoOffDiag(0, 1, 4)
        run.plotdRhoOffDiag(0, 1, 4)
        run.plotRhoFin(0)
        run.plotRhoX(0, 0.1)
        run.plotRhoDiagY(0, 2.5, "-")
        run.plotdRhoDiagY(0, 2.5, "-")
        run.plotRhoOffDiagY(0, 1, 2.5)
        run.plotdRhoOffDiagY(0, 1, 2.5)
        run.plotNeff()
        run.plotEnergyDensity()
        run.plotEntropy()
        run.plotNumberDensity()
        run.plotPArthENoPE()

    def test_example(self):
        """test an example with FortEPiaNORun from explanatory.ini"""
        folder = "output/"
        try:
            os.remove("%s/parthenope.dat" % folder)
            os.remove("%s/parthenope_yi.dat" % folder)
            os.remove("%s/parthenope_rhoee.dat" % folder)
        except:
            pass
        run = fpom.FortEPiaNORun(folder, label="label")
        with open("%s/ini.log" % folder) as _ini:
            ini = _ini.read()
        self.assertEqual(run.ini, ini.replace("\n", " "))
        if "Trh" in ini:
            self.assertEqual(run.zCol, 2)
            self.assertEqual(run.Trhini, 25.0)
            self.assertTrue(run.lowReheating)
        else:
            self.assertEqual(run.zCol, 1)
            self.assertEqual(run.Trhini, None)
            self.assertFalse(run.lowReheating)
        fc = np.loadtxt("%s/fd.dat" % folder)
        self.assertEqualArray(run.yv, fc[:, 0])
        self.assertEqualArray(run.fd, fc[:, 1])
        self.assertEqual(fc.shape[1], 2)
        fc = np.loadtxt("%s/z.dat" % folder)
        Nx = len(fc)
        self.assertEqualArray(run.zdat, fc)
        self.assertEqual(fc.shape[1], 4 if run.lowReheating else 3)
        fc = np.loadtxt("%s/Neff.dat" % folder)
        self.assertEqualArray(run.Neffdat, fc)
        self.assertEqual(fc.shape[1], 4 if run.lowReheating else 3)
        fc = np.loadtxt("%s/energyDensity.dat" % folder)
        self.assertEqualArray(run.endens, fc)
        self.assertEqual(fc.shape[1], 10 if run.lowReheating else 8)
        fc = np.loadtxt("%s/entropy.dat" % folder)
        self.assertEqualArray(run.entropy, fc)
        self.assertEqual(fc.shape[1], 9 if run.lowReheating else 8)
        fc = np.loadtxt("%s/numberDensity.dat" % folder)
        self.assertEqualArray(run.number, fc)
        self.assertEqual(fc.shape[1], 9 if run.lowReheating else 8)
        self.assertTrue(hasattr(run, "resume"))
        self.assertTrue(run.hasResume)
        self.assertIsInstance(run.deltaNeffi, np.ndarray)
        if run.lowReheating:
            self.assertTrue(np.isclose(run.x[0], 0.01, rtol=1e-4))
            self.assertTrue(np.isclose(run.x[-1], 35.0, rtol=1e-4))
            self.assertTrue(np.isclose(run.z[0], 0.01, rtol=1e-4))
            self.assertTrue(np.isclose(run.t[0], 0.001, rtol=1e-4))
            self.assertTrue(np.isclose(run.t[-1], 7.04632596e03, rtol=1e-2))
            self.assertTrue(np.isclose(run.Trh, 25.0, atol=1e-4))
            self.assertTrue(np.isclose(run.Neff, 3.045, atol=1e-3))
            self.assertTrue(np.isclose(run.wfin, 0.6713, atol=1e-3))
            self.assertTrue(np.isclose(run.zfin, 0.9402, atol=1e-3))
            self.assertTrue(np.isclose(run.deltaNeffi[0], 1.016, atol=1e-3))
            self.assertTrue(np.isclose(run.deltaNeffi[1], 1.014, atol=1e-3))
            self.assertTrue(np.isclose(run.deltaNeffi[2], 1.014, atol=1e-3))
        else:
            self.assertTrue(np.isclose(run.x[0], 0.01, rtol=1e-4))
            self.assertTrue(np.isclose(run.x[-1], 35.0, rtol=1e-4))
            self.assertTrue(np.isclose(run.z[0], 1.0288, rtol=1e-4))
            self.assertTrue(np.isclose(run.Neff, 3.0429, atol=1e-4))
            self.assertTrue(np.isclose(run.wfin, 1.0965, atol=1e-4))
            self.assertTrue(np.isclose(run.zfin, 1.5357, atol=1e-4))
            self.assertTrue(np.isclose(run.deltaNeffi[0], 1.0158, atol=1e-4))
            self.assertTrue(np.isclose(run.deltaNeffi[1], 1.0139, atol=1e-4))
            self.assertTrue(np.isclose(run.deltaNeffi[2], 1.0132, atol=1e-4))
        self.assertEqual(len(run.rho), 3)
        self.assertEqual(len(run.rho[0]), 3)
        self.assertEqual(len(run.rho), 3)
        self.assertEqual(len(run.rho[0]), 3)
        for i in range(run.nnu):
            fc = np.loadtxt("%s/nuDens_diag%d.dat" % (folder, i + 1))
            self.assertEqualArray(run.rho[i, i, 0], fc)
            self.assertEqual(run.rho[i, i, 1], None)
            fc = np.loadtxt("%s/nuDens_mass%d.dat" % (folder, i + 1))
            self.assertEqualArray(run.rhoM[i, i, 0], fc)
            self.assertEqual(run.rho[i, i, 1], None)
            for j in range(i + 1, run.nnu):
                fc = np.loadtxt("%s/nuDens_nd_%d%d_re.dat" % (folder, i + 1, j + 1))
                self.assertEqualArray(run.rho[i, j, 0], fc)
                fc = np.loadtxt("%s/nuDens_nd_%d%d_im.dat" % (folder, i + 1, j + 1))
                self.assertEqualArray(run.rho[i, j, 1], fc)
        fc = np.loadtxt("%s/rho_final.dat" % folder)
        self.assertEqualArray(fc.shape, (len(run.yv), 1 + run.nnu))
        fc = np.loadtxt("%s/rho_final_mass.dat" % folder)
        self.assertEqualArray(fc.shape, (len(run.yv), 1 + run.nnu))
        # testing output for BBN
        bbn = np.loadtxt("%s/BBN.dat" % folder)
        self.assertEqualArray(bbn.shape, (Nx, 4))
        rhos = np.loadtxt("%s/rho_tot.dat" % folder)
        self.assertEqualArray(rhos.shape, (Nx, 3 if run.lowReheating else 2))
        fc = np.loadtxt("%s/nuDens_diag1_BBN.dat" % folder)
        self.assertEqualArray(fc.shape, (Nx, 2 + len(run.yv)))
        self.assertTrue(run.hasBBN)
        self.assertEqualArray(run.bbn.shape, (Nx, 4))
        self.assertEqualArray(run.summedrhos.shape, (Nx, 3 if run.lowReheating else 2))
        x99 = next(
            x
            for x, t, r in zip(
                bbn[:, 0], rhos[:, 2 if run.lowReheating else 1], rhos[:, 1]
            )
            if r >= 0.99 * t
        )
        f99 = np.where(bbn[:, 0] >= x99)
        self.assertEqualArray(run.filter99, f99)
        Nx99 = len(run.filter99[0])
        if run.lowReheating:
            self.assertGreater(Nx, Nx99)
        else:
            self.assertEqual(Nx99, Nx)
        self.assertEqualArray(run.parthenope.shape, (Nx99, 7))
        self.assertEqual(
            run.parthenope_cols,
            [
                "x",
                "z",
                "w",
                "rhobarnu",
                "drhobarnu_dx",
                "rho_rad",
                "rho_tot",
            ],
        )
        for f, s in [
            ["parthenope", (Nx99, 7)],
            ["parthenope_yi", (len(run.yv),)],
            ["parthenope_rhoee", (Nx99, len(run.yv))],
        ]:
            fc = np.loadtxt("%s/%s.dat" % (folder, f))
            self.assertEqualArray(fc.shape, s)
        # test that "deltas" creates the delta_ed and delta_nd attributes
        self.assertFalse(hasattr(run, "delta_ed"))
        self.assertFalse(hasattr(run, "delta_nd"))
        run = fpom.FortEPiaNORun(folder, label="label", deltas=True)
        rx = 1 if run.lowReheating else 0
        self.assertEqualArray(
            run.delta_ed,
            [
                (run.endens[-1, run.zCol + 4 + i] - run.endens[rx, run.zCol + 4 + i])
                / run.endens[rx, run.zCol + 4 + i]
                * 100
                for i in range(run.nnu)
            ],
        )
        self.assertEqual(
            run.tot_delta_ed,
            (
                (
                    np.sum(run.endens[-1, run.zCol + 4 : run.zCol + 4 + run.nnu])
                    - np.sum(run.endens[rx, run.zCol + 4 : run.zCol + 4 + run.nnu])
                )
                / np.sum(run.endens[rx, run.zCol + 4 : run.zCol + 4 + run.nnu])
                * 100
            ),
        )
        self.assertEqualArray(
            run.delta_nd,
            [
                (run.number[-1, run.zCol + 4 + i] - run.number[rx, run.zCol + 4 + i])
                / run.number[rx, run.zCol + 4 + i]
                * 100
                for i in range(run.nnu)
            ],
        )
        self.assertEqual(
            run.tot_delta_nd,
            (
                (
                    np.sum(run.number[-1, run.zCol + 4 : run.zCol + 4 + run.nnu])
                    - np.sum(run.number[rx, run.zCol + 4 : run.zCol + 4 + run.nnu])
                )
                / np.sum(run.number[rx, run.zCol + 4 : run.zCol + 4 + run.nnu])
                * 100
            ),
        )
        # now just do plots in order to see that everything works till the end
        self.runAllPlots(run)

    def test_failing(self):
        """test few failing examples with FortEPiaNORun"""
        run = fpom.FortEPiaNORun("output/nonexistent/folder/")
        self.assertFalse(run.lowReheating)
        self.assertFalse(run.hasBBN)
        self.assertFalse(hasattr(run, "yv"))
        self.assertFalse(hasattr(run, "fd"))
        self.assertFalse(hasattr(run, "zdat"))
        self.assertFalse(hasattr(run, "Neffdat"))
        self.assertFalse(hasattr(run, "endens"))
        self.assertFalse(hasattr(run, "entropy"))
        self.assertFalse(hasattr(run, "number"))
        self.assertFalse(hasattr(run, "Neff"))
        self.assertFalse(hasattr(run, "wfin"))
        self.assertFalse(hasattr(run, "zfin"))
        self.assertFalse(hasattr(run, "rho"))
        self.assertFalse(hasattr(run, "rhoM"))
        self.assertFalse(hasattr(run, "resume"))
        self.assertFalse(hasattr(run, "hasResume"))
        self.assertFalse(hasattr(run, "deltaNeffi"))
        self.assertFalse(hasattr(run, "ini"))
        self.assertFalse(hasattr(run, "Trhini"))
        self.assertFalse(hasattr(run, "zCol"))
        self.assertFalse(hasattr(run, "bbn"))
        self.assertFalse(hasattr(run, "parthenope"))
        self.assertFalse(hasattr(run, "filter99"))
        self.runAllPlots(run)

        # repeat creating some bad resume file, e.g. with nans
        folder = "output/no/"
        run = fpom.FortEPiaNORun(folder)
        self.assertTrue(np.isnan(run.yv))
        self.assertTrue(np.isnan(run.fd))
        self.assertEqualArray(run.zdat, np.asarray([[np.nan, np.nan, np.nan]]))
        self.assertEqualArray(run.Neffdat, np.asarray([[np.nan, np.nan, np.nan]]))
        self.assertTrue(np.isnan(run.Neff))
        self.assertTrue(np.isnan(run.wfin))
        self.assertEqual(run.zfin, 1.5)
        self.assertTrue(np.isnan(run.endens))
        self.assertTrue(np.isnan(run.entropy))
        self.assertTrue(np.isnan(run.number))
        self.assertEqual(len(run.rho), 3)
        self.assertEqual(len(run.rho[0]), 3)
        self.assertEqual(len(run.rho), 3)
        self.assertEqual(len(run.rho[0]), 3)
        for i in range(run.nnu):
            self.assertTrue(np.isnan(run.rho[i, i, 0]))
            self.assertTrue(np.isnan(run.rhoM[i, i, 0]))
            for j in range(i + 1, run.nnu):
                self.assertTrue(np.isnan(run.rho[i, j, 0]))
                self.assertTrue(np.isnan(run.rho[i, j, 1]))
        with open("%s/resume.dat" % folder) as _f:
            resume = _f.read()
        self.assertEqual(run.resume, resume.replace("\n", " "))
        self.assertTrue(run.hasResume)
        self.assertEqualArray(run.deltaNeffi, [np.nan, np.nan, np.nan])
        self.assertFalse(run.lowReheating)
        self.assertFalse(run.hasBBN)
        self.runAllPlots(run)

    def test_init(self):
        """test __init__"""
        folder = "output/"
        with patch("fortepianoOutput.FortEPiaNORun.doAllPlots") as _pl, patch(
            "fortepianoOutput.FortEPiaNORun.readIni"
        ) as _ri, patch("fortepianoOutput.FortEPiaNORun.readFD") as _rf, patch(
            "fortepianoOutput.FortEPiaNORun.readWZ"
        ) as _rw, patch(
            "fortepianoOutput.FortEPiaNORun.readNeff"
        ) as _rn, patch(
            "fortepianoOutput.FortEPiaNORun.readEENDensities"
        ) as _rd, patch(
            "fortepianoOutput.FortEPiaNORun.readResume"
        ) as _rr, patch(
            "fortepianoOutput.FortEPiaNORun.readNuDensMatrix"
        ) as _rm, patch(
            "fortepianoOutput.FortEPiaNORun.prepareBBN"
        ) as _bbn, patch(
            "fortepianoOutput.FortEPiaNORun.prepareRhoFinal"
        ) as _prf, patch(
            "fortepianoOutput.FortEPiaNORun.printTableLine"
        ) as _ptl:
            run = fpom.FortEPiaNORun(folder)
            self.assertEqual(_pl.call_count, 0)
            _ri.assert_called_once_with()
            _rf.assert_called_once_with()
            _rw.assert_called_once_with()
            _rn.assert_called_once_with()
            _rd.assert_called_once_with(False)
            _rr.assert_called_once_with()
            _rm.assert_called_once_with(True)
            _bbn.assert_called_once_with()
            _prf.assert_called_once_with()
            _ptl.assert_called_once_with()
        self.assertEqual(run.folder, folder)
        self.assertTrue(run.full)
        self.assertEqual(run.label, "")
        self.assertEqual(run.nnu, 3)
        self.assertTrue(run.verbose)

        with patch("fortepianoOutput.FortEPiaNORun.doAllPlots") as _pl, patch(
            "fortepianoOutput.FortEPiaNORun.readIni"
        ) as _ri, patch("fortepianoOutput.FortEPiaNORun.readFD") as _rf, patch(
            "fortepianoOutput.FortEPiaNORun.readWZ"
        ) as _rw, patch(
            "fortepianoOutput.FortEPiaNORun.readNeff"
        ) as _rn, patch(
            "fortepianoOutput.FortEPiaNORun.readEENDensities"
        ) as _rd, patch(
            "fortepianoOutput.FortEPiaNORun.readResume"
        ) as _rr, patch(
            "fortepianoOutput.FortEPiaNORun.readNuDensMatrix"
        ) as _rm, patch(
            "fortepianoOutput.FortEPiaNORun.prepareBBN"
        ) as _bbn, patch(
            "fortepianoOutput.FortEPiaNORun.prepareRhoFinal"
        ) as _prf, patch(
            "fortepianoOutput.FortEPiaNORun.printTableLine"
        ) as _ptl:
            run = fpom.FortEPiaNORun(
                folder,
                nnu=2,
                label="l",
                deltas=True,
                full=False,
                verbose=False,
                plots=True,
            )
            _pl.assert_called_once_with()
            _ri.assert_called_once_with()
            _rf.assert_called_once_with()
            _rw.assert_called_once_with()
            _rn.assert_called_once_with()
            _rd.assert_called_once_with(True)
            _rr.assert_called_once_with()
            _rm.assert_called_once_with(False)
            _bbn.assert_called_once_with()
            _prf.assert_called_once_with()
            _ptl.assert_called_once_with()
        self.assertEqual(run.folder, folder)
        self.assertFalse(run.full)
        self.assertEqual(run.label, "l")
        self.assertEqual(run.nnu, 2)
        self.assertFalse(run.verbose)

    def test_checkZdat(self):
        """test the checkZdat property"""
        run = fpom.FortEPiaNORun("output/nonexistent")
        self.assertFalse(run.checkZdat())
        run.zdat = "abc"
        self.assertFalse(run.checkZdat())
        run.zdat = np.array([[np.nan, np.nan, np.nan]])
        self.assertFalse(run.checkZdat())
        run.zdat = np.array([[np.nan, np.nan, np.nan], [np.nan, np.nan, np.nan]])
        self.assertFalse(run.checkZdat())
        run.zdat = np.array([[1, 2, 3]])
        self.assertTrue(run.checkZdat())
        run.zdat = np.array([[1, 2, 3], [4, 5, 6]])
        self.assertTrue(run.checkZdat())

    def test_x(self):
        """test the x property"""
        run = fpom.FortEPiaNORun("output/nonexistent")
        with self.assertRaises(AttributeError):
            run.x
        run.zdat = np.array(
            [
                [11.0, 12.0, 13.0, 14.0],
                [12.3, 13.4, 14.5, 15.6],
                [13.4, 14.5, 15.6, 16.7],
            ]
        )
        self.assertEqualArray(run.x, run.zdat[:, 0])
        with patch("fortepianoOutput.FortEPiaNORun.checkZdat", return_value=False):
            with self.assertRaises(AttributeError):
                run.x
            run.bbn = np.array(
                [
                    [1.0, 2.0, 3.0, 4.0],
                    [2.3, 3.4, 4.5, 5.6],
                    [3.4, 4.5, 5.6, 6.7],
                ]
            )
            self.assertEqualArray(run.x, run.bbn[:, 0])

    def test_z(self):
        """test the z property"""
        run = fpom.FortEPiaNORun("output/nonexistent")
        with self.assertRaises(AttributeError):
            run.z
        run.lowReheating = False
        run.zdat = np.array(
            [
                [11.0, 12.0, 13.0, 14.0],
                [12.3, 13.4, 14.5, 15.6],
                [13.4, 14.5, 15.6, 16.7],
            ]
        )
        self.assertEqualArray(run.z, run.zdat[:, 1])
        with patch("fortepianoOutput.FortEPiaNORun.checkZdat", return_value=False):
            with self.assertRaises(AttributeError):
                run.x
            run.bbn = np.array(
                [
                    [1.0, 2.0, 3.0, 4.0],
                    [2.3, 3.4, 4.5, 5.6],
                    [3.4, 4.5, 5.6, 6.7],
                ]
            )
            self.assertEqualArray(run.z, run.bbn[:, 1])
        run.lowReheating = True
        self.assertEqualArray(run.z, run.zdat[:, 2])
        with patch("fortepianoOutput.FortEPiaNORun.checkZdat", return_value=False):
            self.assertEqualArray(run.z, run.bbn[:, 1])

    def test_Tgamma(self):
        """test the Tgamma property"""
        run = fpom.FortEPiaNORun("output/nonexistent")
        run.lowReheating = False
        run.zdat = np.array(
            [
                [11.0, 12.0, 13.0, 14.0],
                [12.3, 13.4, 14.5, 15.6],
                [13.4, 14.5, 15.6, 16.7],
            ]
        )
        self.assertEqualArray(
            run.Tgamma, run.zdat[:, 1] * fpom.ELECTRONMASS_MEV / run.zdat[:, 0]
        )
        run.lowReheating = True
        self.assertEqualArray(
            run.Tgamma, run.zdat[:, 2] * fpom.ELECTRONMASS_MEV / run.zdat[:, 0]
        )

    def test_t(self):
        """test the t property"""
        run = fpom.FortEPiaNORun("output/nonexistent")
        run.lowReheating = False
        run.zdat = np.array(
            [
                [11.0, 12.0, 13.0, 14.0],
                [12.3, 13.4, 14.5, 15.6],
                [13.4, 14.5, 15.6, 16.7],
            ]
        )
        self.assertEqualArray(run.t, np.nan)
        run.lowReheating = True
        self.assertEqualArray(run.t, run.zdat[:, 1])

    def test_drhonu_dx(self):
        """test the drhonu_dx property"""
        run = fpom.FortEPiaNORun("output/nonexistent")
        run.hasBBN = False
        self.assertTrue(np.isnan(run.drhonu_dx))
        run.hasBBN = True
        with self.assertRaises(AttributeError):
            run.drhonu_dx
        run.bbn = np.array(
            [
                [11.0, 12.0, 13.0, 14.0],
                [12.3, 13.4, 14.5, 15.6],
                [13.4, 14.5, 15.6, 16.7],
            ]
        )
        run.zdat = run.bbn
        self.assertEqualArray(run.drhonu_dx, np.gradient(run.bbn[:, 3], run.bbn[:, 0]))

    def test_drhonu_dx_savgol(self):
        """test the drhonu_dx_savgol property"""
        run = fpom.FortEPiaNORun("output/nonexistent")
        run.hasBBN = False
        self.assertTrue(np.isnan(run.drhonu_dx_savgol))
        run.hasBBN = True
        with self.assertRaises(AttributeError):
            run.drhonu_dx_savgol
        run.bbn = np.array(self.explanatory.bbn)
        run.zdat = run.bbn
        self.assertEqualArray(
            run.drhonu_dx_savgol,
            fpom.savgol_filter(np.clip(run.drhonu_dx, 1e-10, None), 51, 1),
        )

    def test_N_func(self):
        """test the N_func property"""
        run = fpom.FortEPiaNORun("output/nonexistent")
        run.hasBBN = False
        self.assertTrue(np.isnan(run.N_func))
        run.hasBBN = True
        with self.assertRaises(AttributeError):
            run.N_func
        run.bbn = np.array(
            [
                [11.0, 12.0, 13.0, 14.0],
                [12.3, 13.4, 14.5, 15.6],
                [13.4, 14.5, 15.6, 16.7],
            ]
        )
        run.zdat = run.bbn
        self.assertEqualArray(
            run.N_func, run.bbn[:, 0] * run.drhonu_dx / run.bbn[:, 1] ** 4
        )

    def test_N_savgol(self):
        """test the N_savgol property"""
        run = fpom.FortEPiaNORun("output/nonexistent")
        run.hasBBN = False
        self.assertTrue(np.isnan(run.N_savgol))
        run.hasBBN = True
        with self.assertRaises(AttributeError):
            run.N_savgol
        run.bbn = np.array(self.explanatory.bbn)
        run.zdat = run.bbn
        self.assertEqualArray(
            run.N_savgol, fpom.savgol_filter(np.clip(run.N_func, 1e-11, None), 75, 1)
        )

    def test_readEENDensities(self):
        """test readEENDensities"""
        run = fpom.FortEPiaNORun("output/nonexistent")
        run.folder = "output/nonexistent1"
        self.assertFalse(hasattr(run, "endens"))
        self.assertFalse(hasattr(run, "entropy"))
        self.assertFalse(hasattr(run, "number"))
        self.assertFalse(hasattr(run, "delta_ed"))
        self.assertFalse(hasattr(run, "tot_delta_ed"))
        self.assertFalse(hasattr(run, "delta_nd"))
        self.assertFalse(hasattr(run, "tot_delta_nd"))
        if os.path.exists(run.folder):
            shutil.rmtree(run.folder)
        os.mkdir(run.folder)
        run.readEENDensities()
        self.assertTrue(np.isnan(run.endens))
        self.assertTrue(np.isnan(run.entropy))
        self.assertTrue(np.isnan(run.number))
        with patch("numpy.loadtxt", return_value=np.array([0, 0])) as _l:
            run.readEENDensities(deltas=True)
            _l.assert_any_call("%s/energyDensity.dat" % run.folder)
            _l.assert_any_call("%s/entropy.dat" % run.folder)
            _l.assert_any_call("%s/numberDensity.dat" % run.folder)
        self.assertEqualArray(run.endens, [0, 0])
        self.assertEqualArray(run.entropy, [0, 0])
        self.assertEqualArray(run.number, [0, 0])
        self.assertFalse(hasattr(run, "delta_ed"))
        self.assertFalse(hasattr(run, "tot_delta_ed"))
        self.assertFalse(hasattr(run, "delta_nd"))
        self.assertFalse(hasattr(run, "tot_delta_nd"))
        run.folder = "output/"
        run.readIni()
        run.readEENDensities()
        self.assertEqualArray(
            run.endens, np.loadtxt("%s/energyDensity.dat" % run.folder)
        )
        self.assertEqualArray(run.entropy, np.loadtxt("%s/entropy.dat" % run.folder))
        self.assertEqualArray(
            run.number, np.loadtxt("%s/numberDensity.dat" % run.folder)
        )
        self.assertFalse(hasattr(run, "delta_ed"))
        self.assertFalse(hasattr(run, "tot_delta_ed"))
        self.assertFalse(hasattr(run, "delta_nd"))
        self.assertFalse(hasattr(run, "tot_delta_nd"))
        run.readEENDensities(deltas=True)
        self.assertEqualArray(
            run.endens, np.loadtxt("%s/energyDensity.dat" % run.folder)
        )
        self.assertEqualArray(run.entropy, np.loadtxt("%s/entropy.dat" % run.folder))
        self.assertEqualArray(
            run.number, np.loadtxt("%s/numberDensity.dat" % run.folder)
        )
        rx = 1 if run.lowReheating else 0
        self.assertEqualArray(
            run.delta_ed,
            [
                (run.endens[-1, run.zCol + 4 + i] - run.endens[rx, run.zCol + 4 + i])
                / run.endens[rx, run.zCol + 4 + i]
                * 100
                for i in range(run.nnu)
            ],
        )
        self.assertEqual(
            run.tot_delta_ed,
            (
                (
                    np.sum(run.endens[-1, run.zCol + 4 : run.zCol + 4 + run.nnu])
                    - np.sum(run.endens[rx, run.zCol + 4 : run.zCol + 4 + run.nnu])
                )
                / np.sum(run.endens[rx, run.zCol + 4 : run.zCol + 4 + run.nnu])
                * 100
            ),
        )
        self.assertEqualArray(
            run.delta_nd,
            [
                (run.number[-1, run.zCol + 4 + i] - run.number[rx, run.zCol + 4 + i])
                / run.number[rx, run.zCol + 4 + i]
                * 100
                for i in range(run.nnu)
            ],
        )
        self.assertEqual(
            run.tot_delta_nd,
            (
                (
                    np.sum(run.number[-1, run.zCol + 4 : run.zCol + 4 + run.nnu])
                    - np.sum(run.number[rx, run.zCol + 4 : run.zCol + 4 + run.nnu])
                )
                / np.sum(run.number[rx, run.zCol + 4 : run.zCol + 4 + run.nnu])
                * 100
            ),
        )

    def test_readFD(self):
        """test readFD"""
        run = fpom.FortEPiaNORun("output/nonexistent")
        for e in [IOError, OSError]:
            with patch("numpy.loadtxt", side_effect=e) as _l:
                run.readFD()
                _l.assert_called_once_with("output/nonexistent/fd.dat")
                self.assertTrue(np.isnan(run.fd))
                self.assertTrue(np.isnan(run.yv))
        with patch("numpy.loadtxt", return_value=np.array([[1, 2], [3, 4]])) as _l:
            run.readFD()
            _l.assert_called_once_with("output/nonexistent/fd.dat")
            self.assertEqualArray(run.fd, [2, 4])
            self.assertEqualArray(run.yv, [1, 3])

    def test_readIni(self):
        """test readIni"""
        run = fpom.FortEPiaNORun("output/nonexistent")
        run.folder = "output/nonexistent1"
        if os.path.exists(run.folder):
            shutil.rmtree(run.folder)
        os.mkdir(run.folder)
        self.assertFalse(hasattr(run, "zCol"))
        self.assertFalse(hasattr(run, "Trhini"))
        self.assertFalse(run.lowReheating)
        run.readIni()
        self.assertEqual(run.ini, "")
        self.assertEqual(run.zCol, 1)
        self.assertEqual(run.Trhini, None)
        self.assertFalse(run.lowReheating)
        with open("%s/ini.log" % run.folder, "w") as _i:
            _i.write("abcd")
        run.readIni()
        self.assertEqual(run.ini, "abcd")
        self.assertEqual(run.zCol, 1)
        self.assertEqual(run.Trhini, None)
        self.assertFalse(run.lowReheating)
        with open("%s/ini.log" % run.folder, "w") as _i:
            _i.write("abcTrh = -1.2E+01def")
        run.readIni()
        self.assertEqual(run.ini, "abcTrh = -1.2E+01def")
        self.assertEqual(run.zCol, 1)
        self.assertEqual(run.Trhini, -12.0)
        self.assertFalse(run.lowReheating)
        with open("%s/ini.log" % run.folder, "w") as _i:
            _i.write("abcTrh = 1.2E-1def")
        run.readIni()
        self.assertEqual(run.ini, "abcTrh = 1.2E-1def")
        self.assertEqual(run.zCol, 2)
        self.assertEqual(run.Trhini, 0.12)
        self.assertTrue(run.lowReheating)

    def test_readNeff(self):
        """test readNeff"""
        run = fpom.FortEPiaNORun("output/nonexistent")
        for l in [False, True]:
            run.lowReheating = l
            for e in [IOError, OSError]:
                with patch("numpy.loadtxt", side_effect=e) as _l:
                    run.readNeff()
                    _l.assert_called_once_with("output/nonexistent/Neff.dat")
                    self.assertEqualArray(
                        run.Neffdat,
                        [[np.nan, np.nan, np.nan, np.nan]]
                        if l
                        else [[np.nan, np.nan, np.nan]],
                    )
            with patch("numpy.loadtxt", return_value=np.array([[1, 2], [3, 4]])) as _l:
                run.readNeff()
                _l.assert_called_once_with("output/nonexistent/Neff.dat")
                self.assertEqualArray(run.Neffdat, [[1, 2], [3, 4]])

    def test_readNuDensMatrix(self):
        """test readNuDensMatrix"""
        run = fpom.FortEPiaNORun("output/nonexistent", nnu=4)
        run.readNuDensMatrix()
        for i in range(run.nnu):
            self.assertTrue(np.isnan(run.rho[i, i, 0]))
            self.assertTrue(np.isnan(run.rhoM[i, i, 0]))
            self.assertEqual(run.rho[i, i, 1], None)
            self.assertEqual(run.rhoM[i, i, 1], None)
            for j in range(i + 1, run.nnu):
                self.assertTrue(np.isnan(run.rho[i, j, 0]))
                self.assertTrue(np.isnan(run.rho[i, j, 1]))
                self.assertEqual(run.rhoM[i, j, 0], None)
                self.assertEqual(run.rhoM[i, j, 1], None)
        run.readNuDensMatrix(full=False)
        for i in range(run.nnu):
            self.assertTrue(np.isnan(run.rho[i, i, 0]))
            self.assertTrue(np.isnan(run.rhoM[i, i, 0]))
            self.assertEqual(run.rho[i, i, 1], None)
            self.assertEqual(run.rhoM[i, i, 1], None)
            for j in range(i + 1, run.nnu):
                self.assertEqual(run.rho[i, j, 0], None)
                self.assertEqual(run.rho[i, j, 1], None)
                self.assertEqual(run.rhoM[i, j, 0], None)
                self.assertEqual(run.rhoM[i, j, 1], None)
        run.folder = "output/"
        run.readNuDensMatrix(full=False)
        for i in range(run.nnu):
            if i < 3:
                self.assertEqualArray(
                    run.rho[i, i, 0],
                    np.loadtxt("%s/nuDens_diag%d.dat" % (run.folder, i + 1)),
                )
                self.assertEqualArray(
                    run.rhoM[i, i, 0],
                    np.loadtxt("%s/nuDens_mass%d.dat" % (run.folder, i + 1)),
                )
            else:
                self.assertTrue(np.isnan(run.rho[i, i, 0]))
                self.assertTrue(np.isnan(run.rhoM[i, i, 0]))
            self.assertEqual(run.rho[i, i, 1], None)
            self.assertEqual(run.rhoM[i, i, 1], None)
            for j in range(i + 1, run.nnu):
                self.assertEqual(run.rho[i, j, 0], None)
                self.assertEqual(run.rho[i, j, 1], None)
                self.assertEqual(run.rhoM[i, j, 0], None)
                self.assertEqual(run.rhoM[i, j, 1], None)
        run.readNuDensMatrix(full=True)
        for i in range(run.nnu):
            if i < 3:
                self.assertEqualArray(
                    run.rho[i, i, 0],
                    np.loadtxt("%s/nuDens_diag%d.dat" % (run.folder, i + 1)),
                )
                self.assertEqualArray(
                    run.rhoM[i, i, 0],
                    np.loadtxt("%s/nuDens_mass%d.dat" % (run.folder, i + 1)),
                )
            else:
                self.assertTrue(np.isnan(run.rho[i, i, 0]))
                self.assertTrue(np.isnan(run.rhoM[i, i, 0]))
            self.assertEqual(run.rho[i, i, 1], None)
            self.assertEqual(run.rhoM[i, i, 1], None)
            for j in range(i + 1, run.nnu):
                if i < 3 and j < 3:
                    self.assertEqualArray(
                        run.rho[i, j, 0],
                        np.loadtxt(
                            "%s/nuDens_nd_%d%d_re.dat" % (run.folder, i + 1, j + 1)
                        ),
                    )
                    self.assertEqualArray(
                        run.rho[i, j, 1],
                        np.loadtxt(
                            "%s/nuDens_nd_%d%d_im.dat" % (run.folder, i + 1, j + 1)
                        ),
                    )
                else:
                    self.assertTrue(np.isnan(run.rho[i, j, 0]))
                    self.assertTrue(np.isnan(run.rho[i, j, 1]))
                self.assertEqual(run.rhoM[i, j, 0], None)
                self.assertEqual(run.rhoM[i, j, 1], None)

    def test_readResume(self):
        """test readResume"""
        run = fpom.FortEPiaNORun("output/nonexistent")
        run.folder = "output/nonexistent2"
        if os.path.exists(run.folder):
            shutil.rmtree(run.folder)
        os.mkdir(run.folder)
        run.readResume()
        self.assertEqual(run.resume, "")
        self.assertFalse(run.hasResume)
        self.assertFalse(hasattr(run, "Neff"))
        self.assertFalse(hasattr(run, "zfin"))
        self.assertFalse(hasattr(run, "wfin"))
        self.assertFalse(hasattr(run, "Trh"))
        self.assertFalse(hasattr(run, "deltaNeffi"))
        t = "Neff=abc"
        with open("%s/resume.dat" % run.folder, "w") as _f:
            _f.write(t)
        run.readResume()
        self.assertEqual(run.resume, t)
        self.assertTrue(run.hasResume)
        self.assertTrue(np.isnan(run.Neff))
        self.assertTrue(np.isnan(run.zfin))
        self.assertTrue(np.isnan(run.wfin))
        self.assertFalse(hasattr(run, "Trh"))
        self.assertEqualArray(run.deltaNeffi, [np.nan for i in range(run.nnu)])
        t = "Neff=3.21\nfinal z =  a2bc"
        with open("%s/resume.dat" % run.folder, "w") as _f:
            _f.write(t)
        run.readResume()
        self.assertEqual(run.resume, t.replace("\n", " "))
        self.assertTrue(run.hasResume)
        self.assertEqual(run.Neff, 3.21)
        self.assertTrue(np.isnan(run.zfin))
        self.assertTrue(np.isnan(run.wfin))
        self.assertFalse(hasattr(run, "Trh"))
        self.assertEqualArray(run.deltaNeffi, [np.nan for i in range(run.nnu)])
        t = "Neff=3.21\nfinal z =  12bc\nfinal w=#"
        with open("%s/resume.dat" % run.folder, "w") as _f:
            _f.write(t)
        run.readResume()
        self.assertEqual(run.resume, t.replace("\n", " "))
        self.assertTrue(run.hasResume)
        self.assertEqual(run.Neff, 3.21)
        self.assertEqual(run.zfin, 12)
        self.assertTrue(np.isnan(run.wfin))
        self.assertFalse(hasattr(run, "Trh"))
        self.assertEqualArray(run.deltaNeffi, [np.nan for i in range(run.nnu)])
        run.lowReheating = False
        t = "Neff=3.21\nfinal z =  12bc\nfinal w =-0.12\nTrh=1.44"
        with open("%s/resume.dat" % run.folder, "w") as _f:
            _f.write(t)
        run.readResume()
        self.assertEqual(run.resume, t.replace("\n", " "))
        self.assertTrue(run.hasResume)
        self.assertEqual(run.Neff, 3.21)
        self.assertEqual(run.zfin, 12)
        self.assertEqual(run.wfin, -0.12)
        self.assertFalse(hasattr(run, "Trh"))
        self.assertEqualArray(run.deltaNeffi, [np.nan for i in range(run.nnu)])
        run.lowReheating = True
        t = "Neff=3.21\nfinal z =  12bc\nfinal w =-0.12\nTrh="
        with open("%s/resume.dat" % run.folder, "w") as _f:
            _f.write(t)
        run.readResume()
        self.assertEqual(run.resume, t.replace("\n", " "))
        self.assertTrue(run.hasResume)
        self.assertEqual(run.Neff, 3.21)
        self.assertEqual(run.zfin, 12)
        self.assertEqual(run.wfin, -0.12)
        self.assertTrue(np.isnan(run.Trh))
        self.assertEqualArray(run.deltaNeffi, [np.nan for i in range(run.nnu)])
        run.Trhini = "0"
        t = "Neff=3.21\nfinal z =  12bc\nfinal w =-0.12\nTrh=1.44"
        with open("%s/resume.dat" % run.folder, "w") as _f:
            _f.write(t)
        with self.assertRaises(ValueError):
            run.readResume()
        self.assertEqual(run.resume, t.replace("\n", " "))
        self.assertTrue(run.hasResume)
        self.assertEqual(run.Neff, 3.21)
        self.assertEqual(run.zfin, 12)
        self.assertEqual(run.wfin, -0.12)
        self.assertEqual(run.Trh, 1.44)
        self.assertEqualArray(run.deltaNeffi, [np.nan for i in range(run.nnu)])
        run.Trhini = 1.440000001
        t = (
            "Neff=3.21\nfinal z =  12bc\nfinal w =-0.12\nTrh=1.44"
            + "\ndeltaneff_0 = 1.015\ndeltaNeff_a = 1.014"
        )
        with open("%s/resume.dat" % run.folder, "w") as _f:
            _f.write(t)
        run.readResume()
        self.assertEqual(run.resume, t.replace("\n", " "))
        self.assertTrue(run.hasResume)
        self.assertEqual(run.Neff, 3.21)
        self.assertEqual(run.zfin, 12)
        self.assertEqual(run.wfin, -0.12)
        self.assertEqual(run.Trh, 1.44)
        self.assertEqualArray(run.deltaNeffi, [np.nan for i in range(run.nnu)])
        t = (
            "Neff=3.21\nfinal z =  12bc\nfinal w =-0.12\nTrh=1.44"
            + "\ndeltaNeff_0 = 1.015\ndeltaNeff_a = 1.014\ndeltaNeff_2 = 1.013"
        )
        with open("%s/resume.dat" % run.folder, "w") as _f:
            _f.write(t)
        run.readResume()
        self.assertEqual(run.resume, t.replace("\n", " "))
        self.assertTrue(run.hasResume)
        self.assertEqual(run.Neff, 3.21)
        self.assertEqual(run.zfin, 12)
        self.assertEqual(run.wfin, -0.12)
        self.assertEqual(run.Trh, 1.44)
        self.assertEqualArray(run.deltaNeffi, [1.015, 1.013])

    def test_readWZ(self):
        """test readWZ"""
        run = fpom.FortEPiaNORun("output/nonexistent")
        for l in [False, True]:
            run.lowReheating = l
            for e in [IOError, OSError]:
                with patch("numpy.loadtxt", side_effect=e) as _l:
                    run.readWZ()
                    _l.assert_called_once_with("output/nonexistent/z.dat")
                    self.assertEqualArray(
                        run.zdat,
                        [[np.nan, np.nan, np.nan, np.nan]]
                        if l
                        else [[np.nan, np.nan, np.nan]],
                    )
            with patch("numpy.loadtxt", return_value=np.array([[1, 2], [3, 4]])) as _l:
                run.readWZ()
                _l.assert_called_once_with("output/nonexistent/z.dat")
                self.assertEqualArray(run.zdat, [[1, 2], [3, 4]])

    def test_prepareBBN_failures(self):
        """test the error management in the prepareBBN function"""
        run = fpom.FortEPiaNORun("output/nonexistent")
        run.yv = np.array([0.1 * i for i in range(10)])
        self.assertFalse(run.hasBBN)
        for res in [
            [True, True, True, False],
            [True, True, False, True],
            [True, False, True, True],
            [False, True, True, True],
        ]:
            with patch("os.path.exists", side_effect=res) as _fe:
                run.prepareBBN()
                self.assertFalse(run.hasBBN)
                self.assertFalse(hasattr(run, "parthenope"))
        rungood = self.explanatory
        run.zdat = np.array([[1, 2, 3.0]])
        run.lowReheating = rungood.lowReheating
        run.yv = np.nan
        with patch("os.path.exists", return_value=True) as _fe:
            with patch("numpy.loadtxt", return_value=np.array([0, 1, 2, 3.0])) as _lt:
                run.prepareBBN()
                _lt.assert_any_call("%s/BBN.dat" % run.folder)
                _lt.assert_any_call("%s/y_grid.dat" % run.folder)
                self.assertEqualArray(run.yv, [0, 1, 2, 3.0])
                run.yv = np.array([[1, np.nan, 2], [3, 4, 5.0]])
                _lt.reset_mock()
                run.prepareBBN()
                _lt.assert_any_call("%s/y_grid.dat" % run.folder)
                self.assertEqualArray(run.yv, [0, 1, 2, 3.0])
        run.yv = np.array([0.1 * i for i in range(10)])
        with patch("os.path.exists", return_value=True) as _fe:
            with patch("numpy.loadtxt", return_value=np.array([0, 1, 2, 3.0])) as _lt:
                run.prepareBBN()
                self.assertFalse(run.hasBBN)
                self.assertEqualArray(run.bbn, np.array([0, 1, 2, 3.0]))
                _lt.assert_called_once_with("%s/BBN.dat" % run.folder)
            _fe.assert_any_call("%s/BBN.dat" % run.folder)
            _fe.assert_any_call("%s/nuDens_diag1_BBN.dat" % run.folder)
            _fe.assert_any_call("%s/rho_tot.dat" % run.folder)
            _fe.assert_any_call("%s/y_grid.dat" % run.folder)
            with patch(
                "numpy.loadtxt",
                side_effect=[
                    rungood.bbn,
                    np.array([0.0, 1]),
                    rungood.bbn,
                    np.array([0.0, 1]),
                ],
            ) as _lt:
                with patch(
                    "fortepianoOutput.FortEPiaNORun.checkZdat", return_value=False
                ) as _cz:
                    run.prepareBBN()
                    self.assertFalse(run.hasBBN)
                with patch(
                    "fortepianoOutput.FortEPiaNORun.checkZdat", return_value=True
                ) as _cz, self.assertRaises(AssertionError):
                    run.prepareBBN()
                    self.assertFalse(run.hasBBN)
            run.zdat = rungood.zdat
            with patch(
                "numpy.loadtxt", side_effect=[rungood.bbn, np.array([0.0, 1])]
            ) as _lt:
                run.prepareBBN()
                self.assertEqual(_lt.call_count, 2)
                _lt.assert_any_call("%s/BBN.dat" % run.folder)
                _lt.assert_any_call("%s/rho_tot.dat" % run.folder)
                if run.lowReheating:
                    self.assertFalse(run.hasBBN)
            with patch(
                "numpy.loadtxt", side_effect=[rungood.bbn, np.array([0.0, 1, 2.0])]
            ) as _lt:
                run.prepareBBN()
                self.assertFalse(run.hasBBN)
            with patch(
                "numpy.loadtxt",
                side_effect=[rungood.bbn, rungood.summedrhos, np.array([0.0, 1, 2.0])],
            ) as _lt, patch("numpy.savetxt") as _st:
                run.prepareBBN()
                self.assertEqual(_lt.call_count, 3)
                _lt.assert_any_call("%s/BBN.dat" % run.folder)
                _lt.assert_any_call("%s/rho_tot.dat" % run.folder)
                _lt.assert_any_call("%s/nuDens_diag1_BBN.dat" % run.folder)
                self.assertTrue(run.hasBBN)
                self.assertEqual(_st.call_count, 0)
            with patch(
                "numpy.loadtxt",
                side_effect=[
                    rungood.bbn,
                    rungood.summedrhos,
                    np.array([[0.0, 1, 2.0], [3.0, 4.0, 5]]),
                ],
            ) as _lt, patch("numpy.savetxt") as _st:
                run.prepareBBN()
                self.assertTrue(run.hasBBN)
                self.assertEqual(_st.call_count, 0)
            with patch(
                "numpy.loadtxt",
                side_effect=[rungood.bbn, rungood.summedrhos, rungood.bbn],
            ) as _lt, patch("numpy.savetxt") as _st:
                run.prepareBBN()
                self.assertTrue(run.hasBBN)
                self.assertEqual(_st.call_count, 0)
        for fer in [
            [True, True, True, True, False, True, True],
            [True, True, True, True, True, False, True],
            [True, True, True, True, True, True, False],
        ]:
            with patch("os.path.exists", side_effect=fer) as _fe:
                with patch(
                    "numpy.loadtxt",
                    side_effect=[rungood.bbn, rungood.summedrhos, rungood.bbn],
                ) as _lt, patch("numpy.savetxt") as _st:
                    run.prepareBBN()
                    self.assertTrue(run.hasBBN)
                    self.assertEqual(_st.call_count, 3)
                    for i, n in enumerate(
                        ["parthenope", "parthenope_yi", "parthenope_rhoee"]
                    ):
                        self.assertEqual(
                            _st.call_args_list[i][0][0], "%s/%s.dat" % (run.folder, n)
                        )
                    for i, n in enumerate(
                        [run.parthenope, run.yv, run.bbn[:, 2:][run.filter99]]
                    ):
                        self.assertEqualArray(_st.call_args_list[i][0][1], n)
                    for i in range(3):
                        self.assertEqual(_st.call_args_list[i][1]["fmt"], "%15.7e")
                    self.assertEqual(
                        _st.call_args_list[0][1]["header"],
                        "".join(["%16s" % s for s in run.parthenope_cols])[3:],
                    )

    def test_prepareRhoFinal(self):
        """test prepareRhoFinal"""
        run = fpom.FortEPiaNORun("output/nonexistent")
        self.assertFalse(hasattr(run, "wfin"))
        with patch("os.path.exists") as _p:
            run.prepareRhoFinal()
            self.assertEqual(_p.call_count, 0)
        run.wfin = np.nan
        with patch("os.path.exists") as _p:
            run.prepareRhoFinal()
            self.assertEqual(_p.call_count, 0)
        run.wfin = 1.1
        self.assertFalse(hasattr(run, "zfin"))
        with patch("os.path.exists") as _p:
            run.prepareRhoFinal()
            self.assertEqual(_p.call_count, 0)
        run.zfin = np.nan
        with patch("os.path.exists") as _p:
            run.prepareRhoFinal()
            self.assertEqual(_p.call_count, 0)
        run.zfin = 1.5
        for pex in [[True, True, True, True, True, True], [False, False]]:
            with patch("os.path.exists", side_effect=pex) as _p, patch(
                "numpy.loadtxt"
            ) as _l:
                run.prepareRhoFinal()
                self.assertEqual(_l.call_count, 0)
        run.wfin = 1.1
        data = np.array([[0, 1, 2], [3, 4, 5], [6, 7, 8]])
        zid = (11.0 / 4.0) ** (1.0 / 3.0)
        var = np.column_stack(
            (
                data[:, 0],
                np.array(
                    [
                        data[:, i] * (np.exp(data[:, 0] / run.wfin) + 1)
                        for i in range(1, data.shape[1])
                    ]
                ).T,
            )
        )
        nor = np.column_stack(
            (
                var[:, 0],
                np.array(
                    [
                        interp1d(var[:, 0], var[:, i], fill_value="extrapolate")(
                            var[:, 0] * run.wfin
                        )
                        * (zid / run.zfin * run.wfin) ** 4
                        for i in range(1, data.shape[1])
                    ]
                ).T,
            )
        )
        with patch(
            "os.path.exists", side_effect=[True, True, False, True, True, False]
        ) as _p, patch("numpy.loadtxt", return_value=data) as _l, patch(
            "numpy.savetxt"
        ) as _s:
            run.prepareRhoFinal()
            for i, fm in enumerate(["rho_final", "rho_final_mass"]):
                _p.assert_any_call("%s/%s.dat" % (run.folder, fm))
                _p.assert_any_call("%s/%s_norm.dat" % (run.folder, fm))
                _p.assert_any_call("%s/%s_var.dat" % (run.folder, fm))
                _l.assert_any_call("%s/%s.dat" % (run.folder, fm))
                self.assertEqual(
                    _s.call_args_list[2 * i][0][0], "%s/%s_var.dat" % (run.folder, fm)
                )
                self.assertEqualArray(_s.call_args_list[2 * i][0][1], var)
                self.assertEqual(_s.call_args_list[2 * i][1], {"fmt": "%15.7e"})
                self.assertEqual(
                    _s.call_args_list[2 * i + 1][0][0],
                    "%s/%s_norm.dat" % (run.folder, fm),
                )
                self.assertEqualArray(_s.call_args_list[2 * i + 1][0][1], nor)
                self.assertEqual(_s.call_args_list[2 * i + 1][1], {"fmt": "%15.7e"})

    def test_interpolateRhoIJ(self):
        """test interpolateRhoIJ"""
        run = fpom.FortEPiaNORun("output/nonexistent")
        run.rho = np.asarray([])
        run.rhoM = np.asarray([])
        with self.assertRaises(IndexError):
            run.interpolateRhoIJ(0, 0, 5)
        run.rho = np.asarray(
            [
                [
                    [None, None],
                    [None, None],
                ],
                [
                    [None, None],
                    [None, None],
                ],
            ],
        )
        run.rho[0, 0, 0] = np.array([[0, 10, 11, 13], [1, 21, 22, 24], [3, 30, 31, 33]])
        run.rho[0, 0, 1] = np.array([[0, 0, 1, 2], [2, 10, 11, 12], [4, 20, 21, 22]])
        run.rho[0, 1, 0] = np.array([[0, 1, 2, 3], [3, 4, 5, 6], [5, 7, 8, 9]])
        run.rho[0, 1, 1] = np.array([[0, 5, 6, 7], [1, 4, 5, 6], [5, 0, 1, 2]])
        run.rho[1, 0, 0] = np.array([[0, 30, 31, 32], [1, 30, 31, 32], [3, 31, 32, 33]])
        run.rho[1, 0, 1] = np.array([[0, 20, 21, 22], [2, 21, 22, 23], [4, 21, 22, 23]])
        with self.assertRaises(AttributeError):
            run.interpolateRhoIJ(0, 0, 5)
        with self.assertRaises(IndexError):
            run.interpolateRhoIJ(0, 0, 5, mass=True)
        with self.assertRaises(TypeError):
            run.interpolateRhoIJ(1, 1, 5)
        run.rhoM = np.array(run.rho)
        run.rhoM[1, 1, 0] = np.array(
            [[0, 21, 22, 23], [3, 24, 25, 26], [5, 27, 28, 29]]
        )
        run.rhoM[1, 1, 1] = np.array(
            [[0, 25, 26, 27], [1, 25, 26, 27], [5, 25, 26, 27]]
        )
        run.yv = np.array([x for x in range(3)])
        self.assertEqualArray(run.interpolateRhoIJ(0, 0, 2), [[0, 1, 3], [13, 24, 33]])
        self.assertEqualArray(
            run.interpolateRhoIJ(0, 0, 2, y2=True),
            [[0, 1, 3], [y * 4 for y in [13, 24, 33]]],
        )
        self.assertEqualArray(
            run.interpolateRhoIJ(0, 1, 1.5, ri=1), [[0, 1, 5], [6.5, 5.5, 1.5]]
        )
        self.assertEqualArray(run.interpolateRhoIJ(1, 0, 1), [[0, 3], [31, 32]])
        self.assertEqualArray(
            run.interpolateRhoIJ(1, 0, 1, ri=1), [[0, 2, 4], [21, 22, 22]]
        )
        self.assertEqualArray(
            run.interpolateRhoIJ(1, 1, 1, ri=1, mass=True), [[0, 5], [26, 26]]
        )

    def test_interpolateRhoIJ_x(self):
        """test interpolateRhoIJ_x"""
        run = fpom.FortEPiaNORun("output/nonexistent")
        with self.assertRaises(AttributeError):
            run.interpolateRhoIJ_x(0, 0, 5)
        run.rho = np.asarray([])
        run.rhoM = np.asarray([])
        with self.assertRaises(AttributeError):
            run.interpolateRhoIJ_x(0, 0, 5)
        run.yv = np.array([x for x in range(3)])
        with self.assertRaises(IndexError):
            run.interpolateRhoIJ_x(0, 0, 5)
        run.rho = np.asarray(
            [
                [
                    [None, None],
                    [None, None],
                ],
                [
                    [None, None],
                    [None, None],
                ],
            ],
        )
        run.rho[0, 0, 0] = np.array([[0, 10, 20, 30], [1, 11, 21, 31], [3, 13, 23, 33]])
        run.rho[0, 0, 1] = np.array([[0, 0, 1, 2], [2, 10, 11, 12], [4, 20, 21, 22]])
        run.rho[0, 1, 0] = np.array([[0, 1, 2, 3], [3, 4, 5, 6], [5, 7, 8, 9]])
        run.rho[0, 1, 1] = np.array([[0, 5, 6, 7], [1, 4, 5, 6], [5, 0, 1, 2]])
        run.rho[1, 0, 0] = np.array([[0, 30, 31, 32], [1, 30, 31, 32], [3, 31, 32, 33]])
        run.rho[1, 0, 1] = np.array([[0, 20, 21, 22], [2, 21, 22, 23], [4, 21, 22, 23]])
        with self.assertRaises(IndexError):
            run.interpolateRhoIJ_x(0, 0, 5, mass=True)
        with self.assertRaises(TypeError):
            run.interpolateRhoIJ_x(1, 1, 5)
        run.rhoM = np.array(run.rho)
        run.rhoM[1, 1, 0] = np.array(
            [[0, 21, 22, 23], [3, 24, 25, 26], [5, 27, 28, 29]]
        )
        run.rhoM[1, 1, 1] = np.array(
            [[0, 25, 26, 27], [1, 25, 26, 27], [5, 25, 26, 27]]
        )
        self.assertEqualArray(
            run.interpolateRhoIJ_x(0, 0, 2), [[0, 1, 2], [12, 22, 32]]
        )
        self.assertEqualArray(
            run.interpolateRhoIJ_x(0, 0, 2, y2=True),
            [[0, 1, 2], [0, 22, 4 * 32]],
        )
        self.assertEqualArray(
            run.interpolateRhoIJ_x(0, 1, 2, ri=1), [[0, 1, 2], [3, 4, 5]]
        )
        self.assertEqualArray(
            run.interpolateRhoIJ_x(1, 0, 2), [[0, 1, 2], [30.5, 31.5, 32.5]]
        )
        self.assertEqualArray(
            run.interpolateRhoIJ_x(1, 1, 2, ri=1, mass=True), [[0, 1, 2], [25, 26, 27]]
        )

    def test_printTableLine(self):
        """test printTableLine"""
        self.emptyRun.printTableLine()
        self.explanatory.printTableLine()
        self.failedRun.printTableLine()

    def test_plotFD(self):
        """test plotFD"""
        with patch("matplotlib.pyplot.plot") as _p:
            self.emptyRun.plotFD()
            self.assertEqual(_p.call_count, 0)
        run = self.explanatory
        with patch("matplotlib.pyplot.plot") as _p, patch(
            "matplotlib.pyplot.xscale"
        ) as _xs, patch("matplotlib.pyplot.yscale") as _ys, patch(
            "matplotlib.pyplot.xlabel"
        ) as _xl, patch(
            "matplotlib.pyplot.ylabel"
        ) as _yl:
            run.plotFD()
            self.assertEqual(_p.call_count, 1)
            self.assertEqualArray(_p.call_args[0], [run.yv, run.fd])
            self.assertEqual(
                _p.call_args[1],
                {
                    "label": "label",
                    "ls": "-",
                    "marker": ".",
                    "c": "k",
                },
            )
            _xl.assert_called_once_with("$y$")
            _yl.assert_called_once_with(r"$y^2 f(y)$")
            _xs.assert_called_once_with("log")
            _ys.assert_called_once_with("log")
            self.explanatory.plotFD(ls=":", lc="r", lab="l", rescale=1.1, fac=2.0)
            self.assertEqualArray(
                _p.call_args[0],
                [
                    run.yv,
                    2.0
                    * run.fd
                    * (np.exp(run.yv) + 1.0)
                    / (np.exp(run.yv / 1.1) + 1.0),
                ],
            )
            self.assertEqual(
                _p.call_args[1],
                {
                    "label": "l",
                    "ls": ":",
                    "marker": ".",
                    "c": "r",
                },
            )
            self.failedRun.plotFD()
            self.assertEqualArray(_p.call_args[0], [np.nan, np.nan])

    def test_plotZ(self):
        """test plotZ"""
        with patch("matplotlib.pyplot.plot") as _p:
            self.emptyRun.plotZ()
            self.assertEqual(_p.call_count, 0)
        run = self.explanatory
        with patch("matplotlib.pyplot.plot") as _p, patch(
            "matplotlib.pyplot.xscale"
        ) as _xs, patch("matplotlib.pyplot.xlabel") as _xl, patch(
            "matplotlib.pyplot.ylabel"
        ) as _yl:
            run.plotZ()
            self.assertEqual(_p.call_count, 1)
            self.assertEqualArray(_p.call_args[0], fpom.stripRepeated(run.zdat, 0, 1))
            self.assertEqual(
                _p.call_args[1],
                {
                    "label": "label",
                    "ls": "-",
                    "c": "k",
                },
            )
            _xl.assert_called_once_with("$x$")
            _yl.assert_called_once_with(r"$z$")
            _xs.assert_called_once_with("log")
            self.explanatory.plotZ(ls=":", lc="r", lab="l")
            self.assertEqual(
                _p.call_args[1],
                {
                    "label": "l",
                    "ls": ":",
                    "c": "r",
                },
            )
            self.failedRun.plotZ()
            self.assertEqualArray(_p.call_args[0], [[np.nan, np.nan, np.nan]] * 2)

    def test_plotW(self):
        """test plotW"""
        with patch("matplotlib.pyplot.plot") as _p:
            self.emptyRun.plotW()
            self.assertEqual(_p.call_count, 0)
        run = self.explanatory
        with patch("matplotlib.pyplot.plot") as _p, patch(
            "matplotlib.pyplot.xscale"
        ) as _xs, patch("matplotlib.pyplot.xlabel") as _xl, patch(
            "matplotlib.pyplot.ylabel"
        ) as _yl:
            run.plotW()
            self.assertEqual(_p.call_count, 1)
            self.assertEqualArray(_p.call_args[0], fpom.stripRepeated(run.zdat, 0, 2))
            self.assertEqual(
                _p.call_args[1],
                {
                    "label": "label",
                    "ls": "-",
                    "c": "k",
                },
            )
            _xl.assert_called_once_with("$x$")
            _yl.assert_called_once_with(r"$w$")
            _xs.assert_called_once_with("log")
            self.explanatory.plotW(ls=":", lc="r", lab="l")
            self.assertEqual(
                _p.call_args[1],
                {
                    "label": "l",
                    "ls": ":",
                    "c": "r",
                },
            )
            self.failedRun.plotW()
            self.assertEqualArray(_p.call_args[0], [[np.nan, np.nan, np.nan]] * 2)

    def test_plotZoverW(self):
        """test plotZoverW"""
        with patch("matplotlib.pyplot.plot") as _p:
            self.emptyRun.plotZoverW()
            self.assertEqual(_p.call_count, 0)
        run = self.explanatory
        with patch("matplotlib.pyplot.plot") as _p, patch(
            "matplotlib.pyplot.xscale"
        ) as _xs, patch("matplotlib.pyplot.xlabel") as _xl, patch(
            "matplotlib.pyplot.ylabel"
        ) as _yl:
            run.plotZoverW()
            self.assertEqual(_p.call_count, 1)
            self.assertEqualArray(
                _p.call_args[0],
                fpom.stripRepeated(
                    np.asarray([[x[0], x[1] / x[2]] for x in run.zdat]), 0, 1
                ),
            )
            self.assertEqual(
                _p.call_args[1],
                {
                    "label": "label",
                    "ls": "-",
                    "c": "k",
                },
            )
            _xl.assert_called_once_with("$x$")
            _yl.assert_called_once_with(r"$z/w$")
            _xs.assert_called_once_with("log")
            self.explanatory.plotZoverW(ls=":", lc="r", lab="l")
            self.assertEqual(
                _p.call_args[1],
                {
                    "label": "l",
                    "ls": ":",
                    "c": "r",
                },
            )
            self.failedRun.plotZoverW()
            self.assertEqualArray(_p.call_args[0], [[np.nan, np.nan, np.nan]] * 2)

    def test_plotDeltaZ(self):
        """test plotDeltaZ"""
        run = self.explanatory
        run1 = fpom.FortEPiaNORun("output/nonexistent")
        run1.zdat = np.array(run.zdat)
        run1.zdat[:, 1] = run.zdat[:, 1] * 0.9
        with patch("matplotlib.pyplot.plot") as _p:
            self.emptyRun.plotDeltaZ(run)
            self.assertEqual(_p.call_count, 0)
            run.plotDeltaZ(self.emptyRun)
            self.assertEqual(_p.call_count, 0)
        with patch("matplotlib.pyplot.plot") as _p, patch(
            "matplotlib.pyplot.xscale"
        ) as _xs, patch("matplotlib.pyplot.xlabel") as _xl, patch(
            "matplotlib.pyplot.ylabel"
        ) as _yl:
            run.plotDeltaZ(run)
            self.assertEqual(_p.call_count, 1)
            xv, yv = fpom.stripRepeated(run.zdat, 0, 1)
            self.assertEqualArray(_p.call_args[0], np.asarray([xv, [0.0 for x in xv]]))
            self.assertEqual(
                _p.call_args[1],
                {
                    "label": "label",
                    "ls": "-",
                    "c": "k",
                },
            )
            _xl.assert_called_once_with("$x$")
            _yl.assert_called_once_with(r"$z-z_{\rm ref}$")
            _xs.assert_called_once_with("log")
            run.plotDeltaZ(run1)
            self.assertEqualArray(_p.call_args[0], np.asarray([xv, -0.1 * yv]))
            run1.plotDeltaZ(run)
            self.assertEqualArray(_p.call_args[0], np.asarray([xv, 0.1 * yv]))
            self.explanatory.plotDeltaZ(run, lab="l", ls=":", lc="r")
            self.assertEqual(
                _p.call_args[1],
                {
                    "label": "l",
                    "ls": ":",
                    "c": "r",
                },
            )
            self.failedRun.plotDeltaZ(run)
            self.assertEqualArray(_p.call_args[0], [[np.nan, np.nan, np.nan]] * 2)
            run.plotDeltaZ(self.failedRun)
            self.assertEqualArray(_p.call_args[0], [xv, [np.nan] * len(xv)])

    def test_plotRhoDiag(self):
        """test plotRhoDiag"""
        run = fpom.FortEPiaNORun("output/nonexistent")
        with patch("matplotlib.pyplot.plot") as _plt:
            run.plotRhoDiag(0, 2, "-")
            run.plotRhoDiag(0, 2, "-", mass=True)
            self.assertEqual(_plt.call_count, 0)
        run.rho = np.asarray([])
        run.rhoM = np.asarray([])
        run.rho = np.asarray(
            [
                [
                    [None, None],
                    [None, None],
                ],
                [
                    [None, None],
                    [None, None],
                ],
            ],
        )
        run.rhoM = np.asarray(
            [
                [
                    [None, None],
                    [None, None],
                ],
                [
                    [None, None],
                    [None, None],
                ],
            ],
        )
        run.yv = np.linspace(0.01, 20, 200)
        run.rho[0, 0, 0] = np.array(
            [[0] + list(run.yv), [1] + list(1.0 / (np.exp(run.yv) + 1))]
        )
        run.rho[1, 1, 0] = np.array(
            [[0] + list(run.yv), [1] + list(2.0 / (np.exp(run.yv) + 1))]
        )
        with patch("matplotlib.pyplot.plot") as _plt, patch(
            "matplotlib.pyplot.xscale"
        ) as _xs, patch("matplotlib.pyplot.xlabel") as _xl, patch(
            "matplotlib.pyplot.ylabel"
        ) as _yl:
            run.plotRhoDiag(0, 1, "-")
            _xs.assert_called_once_with("log")
            _xl.assert_called_once_with("$x$")
            _yl.assert_called_once_with(r"$\rho_{\alpha\alpha}$")
            _plt.assert_called_once()
            self.assertEqualArray(
                _plt.call_args[0], fpom.stripRepeated(run.rho[0, 0, 0], 0, 1)
            )
            self.assertEqual(
                _plt.call_args[1],
                {"label": r"%s $\alpha$=%d" % (run.label, 1), "ls": "-", "c": "k"},
            )
        with patch("matplotlib.pyplot.plot") as _plt:
            run.plotRhoDiag(1, 0, ":", lc="r")
            self.assertEqualArray(
                _plt.call_args[0], fpom.stripRepeated(run.rho[1, 1, 0], 0, 0)
            )
            self.assertEqual(
                _plt.call_args[1],
                {"label": r"%s $\alpha$=%d" % (run.label, 2), "ls": ":", "c": "r"},
            )
        run.rhoM[1, 1, 0] = np.array(
            [
                [0] + list(1.0 / (np.exp(run.yv) + 1)),
                [1] + list(2.0 / (np.exp(run.yv) + 1)),
            ]
        )
        with patch("matplotlib.pyplot.plot") as _plt:
            run.plotRhoDiag(1, 1, "-", mass=True)
            self.assertEqualArray(
                _plt.call_args[0], fpom.stripRepeated(run.rhoM[1, 1, 0], 0, 1)
            )
            self.assertEqual(
                _plt.call_args[1],
                {"label": r"%s $\alpha$=%d" % (run.label, 2), "ls": "-", "c": "k"},
            )

    def test_plotdRhoDiag(self):
        """test plotdRhoDiag"""
        run = fpom.FortEPiaNORun("output/nonexistent")
        with patch("matplotlib.pyplot.plot") as _plt:
            run.plotdRhoDiag(0, 2, "-")
            run.plotdRhoDiag(0, 2, "-", mass=True)
            self.assertEqual(_plt.call_count, 0)
        run.rho = np.asarray([])
        run.rhoM = np.asarray([])
        run.rho = np.asarray(
            [
                [
                    [None, None],
                    [None, None],
                ],
                [
                    [None, None],
                    [None, None],
                ],
            ],
        )
        run.rhoM = np.asarray(
            [
                [
                    [None, None],
                    [None, None],
                ],
                [
                    [None, None],
                    [None, None],
                ],
            ],
        )
        run.yv = np.linspace(0.01, 20, 200)
        run.rho[0, 0, 0] = np.array(
            [[0] + list(run.yv), [1] + list(1.0 / (np.exp(run.yv) + 1))]
        )
        run.rho[1, 1, 0] = np.array(
            [[0] + list(run.yv), [1] + list(2.0 / (np.exp(run.yv) + 1))]
        )
        with patch("matplotlib.pyplot.plot") as _plt, patch(
            "matplotlib.pyplot.xscale"
        ) as _xs, patch("matplotlib.pyplot.xlabel") as _xl, patch(
            "matplotlib.pyplot.ylabel"
        ) as _yl:
            run.plotdRhoDiag(0, 1, "-")
            _xs.assert_called_once_with("log")
            _xl.assert_called_once_with("$x$")
            _yl.assert_called_once_with(r"$d\rho_{\alpha\alpha}/dx$")
            _plt.assert_called_once()
            xv, yv = fpom.stripRepeated(run.rho[0, 0, 0], 0, 1)
            self.assertEqualArray(_plt.call_args[0], [xv, np.gradient(yv, xv)])
            self.assertEqual(
                _plt.call_args[1],
                {"label": r"%s $\alpha$=%d" % (run.label, 1), "ls": "-", "c": "k"},
            )
        with patch("matplotlib.pyplot.plot") as _plt:
            run.plotdRhoDiag(1, 0, ":", lc="r")
            xv, yv = fpom.stripRepeated(run.rho[1, 1, 0], 0, 0)
            self.assertEqualArray(_plt.call_args[0], [xv, np.gradient(yv, xv)])
            self.assertEqual(
                _plt.call_args[1],
                {"label": r"%s $\alpha$=%d" % (run.label, 2), "ls": ":", "c": "r"},
            )
        run.rhoM[1, 1, 0] = np.array(
            [
                [0] + list(1.0 / (np.exp(run.yv) + 1)),
                [1] + list(2.0 / (np.exp(run.yv) + 1)),
            ]
        )
        with patch("matplotlib.pyplot.plot") as _plt:
            run.plotdRhoDiag(1, 1, "-", mass=True)
            xv, yv = fpom.stripRepeated(run.rhoM[1, 1, 0], 0, 1)
            self.assertEqualArray(_plt.call_args[0], [xv, np.gradient(yv, xv)])
            self.assertEqual(
                _plt.call_args[1],
                {"label": r"%s $\alpha$=%d" % (run.label, 2), "ls": "-", "c": "k"},
            )

    def test_plotRhoOffDiag(self):
        """test plotRhoOffDiag"""
        run = fpom.FortEPiaNORun("output/nonexistent")
        run.full = True
        with patch("matplotlib.pyplot.plot") as _plt:
            run.plotRhoOffDiag(0, 1, 1)
            run.plotRhoOffDiag(0, 1, 1, mass=True)
            self.assertEqual(_plt.call_count, 0)
        run.rho = np.asarray([])
        run.rhoM = np.asarray([])
        run.rho = np.asarray(
            [
                [
                    [None, None],
                    [None, None],
                ],
                [
                    [None, None],
                    [None, None],
                ],
            ],
        )
        run.rhoM = np.asarray(
            [
                [
                    [None, None],
                    [None, None],
                ],
                [
                    [None, None],
                    [None, None],
                ],
            ],
        )
        run.yv = np.linspace(0.01, 20, 200)
        run.rho[0, 1, 0] = np.array(
            [[0] + list(run.yv), [1] + list(1.0 / (np.exp(run.yv) + 1))]
        )
        run.rho[0, 1, 1] = np.array(
            [[0] + list(run.yv), [1] + list(2.0 / (np.exp(run.yv) + 1))]
        )
        del run.full
        with patch("matplotlib.pyplot.plot") as _plt:
            run.plotRhoOffDiag(0, 1, 1)
            run.plotRhoOffDiag(0, 1, 1, mass=True)
            self.assertEqual(_plt.call_count, 0)
        run.full = True
        with patch("matplotlib.pyplot.plot") as _plt, patch(
            "matplotlib.pyplot.xscale"
        ) as _xs, patch("matplotlib.pyplot.xlabel") as _xl, patch(
            "matplotlib.pyplot.ylabel"
        ) as _yl:
            run.plotRhoOffDiag(0, 1, 1, im=False)
            _xs.assert_called_once_with("log")
            _xl.assert_called_once_with("$x$")
            _yl.assert_called_once_with(r"$\rho_{\alpha\beta}$")
            _plt.assert_called_once()
            xv, yv = fpom.stripRepeated(run.rho[0, 1, 0], 0, 1)
            self.assertEqualArray(_plt.call_args[0], [xv, yv])
            self.assertEqual(
                _plt.call_args[1],
                {
                    "label": r"%s $\alpha\beta$=%d%d re" % (run.label, 1, 2),
                    "ls": "-",
                    "c": "k",
                },
            )
            run.plotRhoOffDiag(0, 1, 1, lc="r")
            self.assertEqual(_plt.call_count, 3)
            xv, yv = fpom.stripRepeated(run.rho[0, 1, 1], 0, 1)
            self.assertEqualArray(_plt.call_args[0], [xv, yv])
            self.assertEqual(
                _plt.call_args[1],
                {
                    "label": r"%s $\alpha\beta$=%d%d im" % (run.label, 1, 2),
                    "ls": ":",
                    "c": "r",
                },
            )
        with patch("matplotlib.pyplot.plot") as _plt:
            run.plotRhoOffDiag(0, 1, 1, lc="r", im=False)
            xv, yv = fpom.stripRepeated(run.rho[0, 1, 0], 0, 1)
            self.assertEqualArray(_plt.call_args[0], [xv, yv])
            self.assertEqual(
                _plt.call_args[1],
                {
                    "label": r"%s $\alpha\beta$=%d%d re" % (run.label, 1, 2),
                    "ls": "-",
                    "c": "r",
                },
            )
        run.rhoM[1, 0, 0] = np.array(
            [
                [0] + list(1.0 / (np.exp(run.yv) + 1)),
                [1] + list(2.0 / (np.exp(run.yv) + 1)),
            ]
        )
        with patch("matplotlib.pyplot.plot") as _plt:
            run.plotRhoOffDiag(1, 0, 0, im=False, mass=True)
            xv, yv = fpom.stripRepeated(run.rhoM[1, 0, 0], 0, 0)
            self.assertEqualArray(_plt.call_args[0], [xv, yv])

    def test_plotdRhoOffDiag(self):
        """test plotdRhoOffDiag"""
        run = fpom.FortEPiaNORun("output/nonexistent")
        run.full = True
        with patch("matplotlib.pyplot.plot") as _plt:
            run.plotdRhoOffDiag(0, 1, 1)
            run.plotdRhoOffDiag(0, 1, 1, mass=True)
            self.assertEqual(_plt.call_count, 0)
        run.rho = np.asarray([])
        run.rhoM = np.asarray([])
        run.rho = np.asarray(
            [
                [
                    [None, None],
                    [None, None],
                ],
                [
                    [None, None],
                    [None, None],
                ],
            ],
        )
        run.rhoM = np.asarray(
            [
                [
                    [None, None],
                    [None, None],
                ],
                [
                    [None, None],
                    [None, None],
                ],
            ],
        )
        run.yv = np.linspace(0.01, 20, 200)
        run.rho[0, 1, 0] = np.array(
            [[0] + list(run.yv), [1] + list(1.0 / (np.exp(run.yv) + 1))]
        )
        run.rho[0, 1, 1] = np.array(
            [[0] + list(run.yv), [1] + list(2.0 / (np.exp(run.yv) + 1))]
        )
        del run.full
        with patch("matplotlib.pyplot.plot") as _plt:
            run.plotdRhoOffDiag(0, 1, 1)
            run.plotdRhoOffDiag(0, 1, 1, mass=True)
            self.assertEqual(_plt.call_count, 0)
        run.full = True
        with patch("matplotlib.pyplot.plot") as _plt, patch(
            "matplotlib.pyplot.xscale"
        ) as _xs, patch("matplotlib.pyplot.xlabel") as _xl, patch(
            "matplotlib.pyplot.ylabel"
        ) as _yl:
            run.plotdRhoOffDiag(0, 1, 1, im=False)
            _xs.assert_called_once_with("log")
            _xl.assert_called_once_with("$x$")
            _yl.assert_called_once_with(r"$d\rho_{\alpha\beta}/dx$")
            _plt.assert_called_once()
            xv, yv = fpom.stripRepeated(run.rho[0, 1, 0], 0, 1)
            self.assertEqualArray(_plt.call_args[0], [xv, np.gradient(yv, xv)])
            self.assertEqual(
                _plt.call_args[1],
                {
                    "label": r"%s $\alpha\beta$=%d%d re" % (run.label, 1, 2),
                    "ls": "-",
                    "c": "k",
                },
            )
            run.plotdRhoOffDiag(0, 1, 1, lc="r")
            self.assertEqual(_plt.call_count, 3)
            xv, yv = fpom.stripRepeated(run.rho[0, 1, 1], 0, 1)
            self.assertEqualArray(_plt.call_args[0], [xv, np.gradient(yv, xv)])
            self.assertEqual(
                _plt.call_args[1],
                {
                    "label": r"%s $\alpha\beta$=%d%d im" % (run.label, 1, 2),
                    "ls": ":",
                    "c": "r",
                },
            )
        with patch("matplotlib.pyplot.plot") as _plt:
            run.plotdRhoOffDiag(0, 1, 1, lc="r", im=False)
            xv, yv = fpom.stripRepeated(run.rho[0, 1, 0], 0, 1)
            self.assertEqualArray(_plt.call_args[0], [xv, np.gradient(yv, xv)])
            self.assertEqual(
                _plt.call_args[1],
                {
                    "label": r"%s $\alpha\beta$=%d%d re" % (run.label, 1, 2),
                    "ls": "-",
                    "c": "r",
                },
            )
        run.rhoM[1, 0, 0] = np.array(
            [
                [0] + list(1.0 / (np.exp(run.yv) + 1)),
                [1] + list(2.0 / (np.exp(run.yv) + 1)),
            ]
        )
        with patch("matplotlib.pyplot.plot") as _plt:
            run.plotdRhoOffDiag(1, 0, 0, im=False, mass=True)
            xv, yv = fpom.stripRepeated(run.rhoM[1, 0, 0], 0, 0)
            self.assertEqualArray(_plt.call_args[0], [xv, np.gradient(yv, xv)])

    def test_plotRhoFin(self):
        """test plotRhoFin"""
        run = fpom.FortEPiaNORun("output/nonexistent")
        with patch("matplotlib.pyplot.plot") as _plt:
            run.plotRhoFin(0)
            run.plotRhoFin(0, mass=True)
            self.assertEqual(_plt.call_count, 0)
        run.rho = np.asarray([])
        run.rhoM = np.asarray([])
        run.rho = np.asarray(
            [
                [
                    [None, None],
                    [None, None],
                ],
                [
                    [None, None],
                    [None, None],
                ],
            ],
        )
        run.rhoM = np.asarray(
            [
                [
                    [None, None],
                    [None, None],
                ],
                [
                    [None, None],
                    [None, None],
                ],
            ],
        )
        run.yv = np.linspace(0.01, 20, 200)
        run.rho[0, 0, 0] = np.array(
            [[0] + list(run.yv), [1] + list(1.0 / (np.exp(run.yv) + 1))]
        )
        run.rho[1, 0, 1] = np.array(
            [[0] + list(run.yv), [1] + list(2.0 / (np.exp(run.yv) + 1))]
        )
        with patch("matplotlib.pyplot.plot") as _plt, patch(
            "matplotlib.pyplot.xscale"
        ) as _xs, patch("matplotlib.pyplot.xlabel") as _xl, patch(
            "matplotlib.pyplot.ylabel"
        ) as _yl:
            run.plotRhoFin(0)
            _xl.assert_called_once_with("$y$")
            _yl.assert_called_once_with(r"$\rho_{\alpha\beta}^{\rm fin}(y)$")
            _plt.assert_called_once()
            self.assertEqualArray(_plt.call_args[0], [run.yv, run.rho[0, 0, 0][-1, 1:]])
            self.assertEqual(
                _plt.call_args[1],
                {
                    "label": r"%s $\alpha\beta$=%d%d %s" % (run.label, 1, 1, "re"),
                    "ls": "-",
                    "c": "k",
                },
            )
        with patch("matplotlib.pyplot.plot") as _plt, patch(
            "matplotlib.pyplot.ylabel"
        ) as _yl:
            run.plotRhoFin(1, i2=0, ri=1, ls=":", lc="r", y2=True)
            _yl.assert_called_once_with(r"$y^2\rho_{\alpha\beta}^{\rm fin}(y)$")
            self.assertEqualArray(
                _plt.call_args[0], [run.yv, run.yv ** 2 * run.rho[1, 0, 1][-1, 1:]]
            )
            self.assertEqual(
                _plt.call_args[1],
                {
                    "label": r"%s $\alpha\beta$=%d%d %s" % (run.label, 2, 1, "im"),
                    "ls": ":",
                    "c": "r",
                },
            )
        run.rhoM[1, 1, 0] = np.array(
            [
                [0] + list(1.0 / (np.exp(run.yv) + 1)),
                [1] + list(2.0 / (np.exp(run.yv) + 1)),
            ]
        )
        with patch("matplotlib.pyplot.plot") as _plt:
            run.plotRhoFin(1, mass=True, lab="lab")
            self.assertEqualArray(
                _plt.call_args[0], [run.yv, run.rhoM[1, 1, 0][-1, 1:]]
            )
            self.assertEqual(_plt.call_args[1], {"label": r"lab", "ls": "-", "c": "k"})

    def test_plotRhoX(self):
        """test plotRhoX"""
        run = fpom.FortEPiaNORun("output/nonexistent")
        with patch("matplotlib.pyplot.plot") as _plt:
            run.plotRhoX(0, 0.5)
            run.plotRhoX(0, 0.5, mass=True)
            self.assertEqual(_plt.call_count, 0)
        run.rho = np.asarray([])
        run.rhoM = np.asarray([])
        run.rho = np.asarray(
            [
                [
                    [None, None],
                    [None, None],
                ],
                [
                    [None, None],
                    [None, None],
                ],
            ],
        )
        run.rhoM = np.asarray(
            [
                [
                    [None, None],
                    [None, None],
                ],
                [
                    [None, None],
                    [None, None],
                ],
            ],
        )
        run.yv = np.linspace(0.01, 20, 200)
        run.rho[0, 0, 0] = np.array(
            [[0] + list(run.yv), [1] + list(1.0 / (np.exp(run.yv) + 1))]
        )
        run.rho[1, 0, 1] = np.array(
            [[0] + list(run.yv), [1] + list(2.0 / (np.exp(run.yv) + 1))]
        )
        with patch("matplotlib.pyplot.plot") as _plt, patch(
            "matplotlib.pyplot.xscale"
        ) as _xs, patch("matplotlib.pyplot.xlabel") as _xl, patch(
            "matplotlib.pyplot.ylabel"
        ) as _yl:
            run.plotRhoX(0, 0.5)
            _xl.assert_called_once_with("$y$")
            _yl.assert_called_once_with(r"$\rho_{\alpha\beta}(y)$")
            _plt.assert_called_once()
            self.assertEqualArray(
                _plt.call_args[0],
                run.interpolateRhoIJ_x(0, 0, 0.5, 0, y2=False, mass=False),
            )
            self.assertEqual(
                _plt.call_args[1],
                {
                    "label": r"%s $\alpha\beta$=%d%d %s x=%f"
                    % (run.label, 1, 1, "re", 0.5),
                    "ls": "-",
                    "c": "k",
                },
            )
        with patch("matplotlib.pyplot.plot") as _plt, patch(
            "matplotlib.pyplot.ylabel"
        ) as _yl:
            run.plotRhoX(1, 0.5, i2=0, ri=1, ls=":", lc="r", y2=True, lab="aaa")
            _yl.assert_called_once_with(r"$y^2\rho_{\alpha\beta}(y)$")
            self.assertEqualArray(
                _plt.call_args[0],
                run.interpolateRhoIJ_x(1, 0, 0.5, 1, y2=True, mass=False),
            )
            self.assertEqual(
                _plt.call_args[1],
                {
                    "label": r"aaa",
                    "ls": ":",
                    "c": "r",
                },
            )
        run.rhoM[1, 1, 0] = np.array(
            [
                [0] + list(1.0 / (np.exp(run.yv) + 1)),
                [1] + list(2.0 / (np.exp(run.yv) + 1)),
            ]
        )
        with patch("matplotlib.pyplot.plot") as _plt:
            run.plotRhoX(1, 0.5, mass=True)
            self.assertEqualArray(
                _plt.call_args[0],
                run.interpolateRhoIJ_x(1, 1, 0.5, 0, y2=False, mass=True),
            )
        with patch("matplotlib.pyplot.plot") as _plt:
            run.plotRhoX(1, 0.5, mass=True, divide_by=2.0)
            x, y = run.interpolateRhoIJ_x(1, 1, 0.5, 0, y2=False, mass=True)
            self.assertEqualArray(
                _plt.call_args[0],
                [x, np.array(y) / 2.0],
            )
        with self.assertRaises(AttributeError):
            run.plotRhoX(1, 0.5, mass=True, divide_fd=True)
        run.fd = np.linspace(0.1, 10, 200)
        with patch("matplotlib.pyplot.plot") as _plt:
            run.plotRhoX(1, 0.5, mass=True, divide_fd=True)
            x, y = run.interpolateRhoIJ_x(1, 1, 0.5, 0, y2=False, mass=True)
            self.assertEqualArray(
                _plt.call_args[0],
                [x, np.array(y) / run.fd],
            )
        with patch("matplotlib.pyplot.plot") as _plt:
            run.plotRhoX(1, 0.5, mass=True, divide_fd=True, divide_by=2.0)
            x, y = run.interpolateRhoIJ_x(1, 1, 0.5, 0, y2=False, mass=True)
            self.assertEqualArray(
                _plt.call_args[0],
                [x, np.array(y) / run.fd / 2.0],
            )

    def test_plotRhoDiagY(self):
        """test plotRhoDiagY"""
        run = fpom.FortEPiaNORun("output/nonexistent")
        with patch("matplotlib.pyplot.plot") as _plt:
            run.plotRhoDiagY(0, 2.5, "-")
            run.plotRhoDiagY(0, 2.5, "-", mass=True)
            self.assertEqual(_plt.call_count, 0)
        run.rho = np.asarray([])
        run.rhoM = np.asarray([])
        run.rho = np.asarray(
            [
                [
                    [None, None],
                    [None, None],
                ],
                [
                    [None, None],
                    [None, None],
                ],
            ],
        )
        run.rhoM = np.asarray(
            [
                [
                    [None, None],
                    [None, None],
                ],
                [
                    [None, None],
                    [None, None],
                ],
            ],
        )
        run.yv = np.linspace(0.01, 20, 200)
        run.rho[0, 0, 0] = np.array(
            [[0] + list(run.yv), [1] + list(1.0 / (np.exp(run.yv) + 1))]
        )
        run.rho[1, 1, 0] = np.array(
            [[0] + list(run.yv), [1] + list(2.0 / (np.exp(run.yv) + 1))]
        )
        with patch("matplotlib.pyplot.plot") as _plt, patch(
            "matplotlib.pyplot.xscale"
        ) as _xs, patch("matplotlib.pyplot.xlabel") as _xl, patch(
            "matplotlib.pyplot.ylabel"
        ) as _yl:
            run.plotRhoDiagY(0, 2.5, "-")
            _xs.assert_called_once_with("log")
            _xl.assert_called_once_with("$x$")
            _yl.assert_called_once_with(r"$\rho_{\alpha\alpha}$")
            _plt.assert_called_once()
            x, yv = run.interpolateRhoIJ(0, 0, 2.5)
            self.assertEqualArray(_plt.call_args[0], [x, yv])
            self.assertEqual(
                _plt.call_args[1],
                {"label": r"%s $\alpha$=%d" % (run.label, 1), "ls": "-", "c": "k"},
            )
        with patch("matplotlib.pyplot.plot") as _plt, patch(
            "matplotlib.pyplot.ylabel"
        ) as _yl:
            run.plotRhoDiagY(1, 2.5, ":", lc="r", y2=True)
            _yl.assert_called_once_with(r"$y^2\rho_{\alpha\alpha}$")
            x, yv = run.interpolateRhoIJ(1, 1, 2.5)
            self.assertEqualArray(_plt.call_args[0], [x, np.array(yv) * 2.5 ** 2])
            self.assertEqual(
                _plt.call_args[1],
                {"label": r"%s $\alpha$=%d" % (run.label, 2), "ls": ":", "c": "r"},
            )
        run.rhoM[1, 1, 0] = np.array(
            [
                [0] + list(1.0 / (np.exp(run.yv) + 1)),
                [1] + list(2.0 / (np.exp(run.yv) + 1)),
            ]
        )
        with patch("matplotlib.pyplot.plot") as _plt:
            run.plotRhoDiagY(1, 2.5, "-", mass=True, lab="lab")
            x, yv = run.interpolateRhoIJ(1, 1, 2.5, mass=True)
            self.assertEqualArray(_plt.call_args[0], [x, yv])
            self.assertEqual(
                _plt.call_args[1],
                {"label": "lab", "ls": "-", "c": "k"},
            )
        with patch("matplotlib.pyplot.plot") as _plt:
            run.plotRhoDiagY(1, 2.5, "-", mass=True, lab="lab", divide_by=3.0)
            x, yv = run.interpolateRhoIJ(1, 1, 2.5, mass=True)
            self.assertEqualArray(_plt.call_args[0], [x, np.asarray(yv) / 3.0])
            self.assertEqual(
                _plt.call_args[1],
                {"label": "lab", "ls": "-", "c": "k"},
            )
        with self.assertRaises(TypeError):
            run.plotRhoDiagY(1, 2.5, "-", mass=True, lab="lab", divide_by="a")

    def test_plotdRhoDiagY(self):
        """test plotdRhoDiagY"""
        run = fpom.FortEPiaNORun("output/nonexistent")
        with patch("matplotlib.pyplot.plot") as _plt:
            run.plotdRhoDiagY(0, 2.5, "-")
            run.plotdRhoDiagY(0, 2.5, "-", mass=True)
            self.assertEqual(_plt.call_count, 0)
        run.rho = np.asarray([])
        run.rhoM = np.asarray([])
        run.rho = np.asarray(
            [
                [
                    [None, None],
                    [None, None],
                ],
                [
                    [None, None],
                    [None, None],
                ],
            ],
        )
        run.rhoM = np.asarray(
            [
                [
                    [None, None],
                    [None, None],
                ],
                [
                    [None, None],
                    [None, None],
                ],
            ],
        )
        run.yv = np.linspace(0.01, 20, 200)
        run.rho[0, 0, 0] = np.array(
            [[0] + list(run.yv), [1] + list(1.0 / (np.exp(run.yv) + 1))]
        )
        run.rho[1, 1, 0] = np.array(
            [[0] + list(run.yv), [1] + list(2.0 / (np.exp(run.yv) + 1))]
        )
        with patch("matplotlib.pyplot.plot") as _plt, patch(
            "matplotlib.pyplot.xscale"
        ) as _xs, patch("matplotlib.pyplot.xlabel") as _xl, patch(
            "matplotlib.pyplot.ylabel"
        ) as _yl:
            run.plotdRhoDiagY(0, 2.5, "-")
            _xs.assert_called_once_with("log")
            _xl.assert_called_once_with("$x$")
            _yl.assert_called_once_with(r"$d\rho_{\alpha\alpha}/dx$")
            _plt.assert_called_once()
            x, yv = run.interpolateRhoIJ(0, 0, 2.5)
            self.assertEqualArray(_plt.call_args[0], [x, np.gradient(yv, x)])
            self.assertEqual(
                _plt.call_args[1],
                {"label": r"%s $\alpha$=%d" % (run.label, 1), "ls": "-", "c": "k"},
            )
        with patch("matplotlib.pyplot.plot") as _plt, patch(
            "matplotlib.pyplot.ylabel"
        ) as _yl:
            run.plotdRhoDiagY(1, 2.5, ":", lc="r", y2=True)
            _yl.assert_called_once_with(r"$dy^2\rho_{\alpha\alpha}/dx$")
            x, yv = run.interpolateRhoIJ(1, 1, 2.5)
            self.assertEqualArray(
                _plt.call_args[0], [x, np.gradient(np.array(yv) * 2.5 ** 2, x)]
            )
            self.assertEqual(
                _plt.call_args[1],
                {"label": r"%s $\alpha$=%d" % (run.label, 2), "ls": ":", "c": "r"},
            )
        run.rhoM[1, 1, 0] = np.array(
            [
                [0] + list(1.0 / (np.exp(run.yv) + 1)),
                [1] + list(2.0 / (np.exp(run.yv) + 1)),
            ]
        )
        with patch("matplotlib.pyplot.plot") as _plt:
            run.plotdRhoDiagY(1, 2.5, "-", mass=True, lab="lab")
            x, yv = run.interpolateRhoIJ(1, 1, 2.5, mass=True)
            self.assertEqualArray(_plt.call_args[0], [x, np.gradient(yv, x)])
            self.assertEqual(
                _plt.call_args[1],
                {"label": "lab", "ls": "-", "c": "k"},
            )

    def test_plotRhoOffDiagY(self):
        """test plotRhoOffDiagY"""
        run = fpom.FortEPiaNORun("output/nonexistent")
        run.full = True
        with patch("matplotlib.pyplot.plot") as _plt:
            run.plotRhoOffDiag(0, 1, 1)
            run.plotRhoOffDiag(0, 1, 1, mass=True)
            self.assertEqual(_plt.call_count, 0)
        run.rho = np.asarray([])
        run.rhoM = np.asarray([])
        run.rho = np.asarray(
            [
                [
                    [None, None],
                    [None, None],
                ],
                [
                    [None, None],
                    [None, None],
                ],
            ],
        )
        run.rhoM = np.asarray(
            [
                [
                    [None, None],
                    [None, None],
                ],
                [
                    [None, None],
                    [None, None],
                ],
            ],
        )
        run.yv = np.linspace(0.01, 20, 200)
        run.rho[0, 1, 0] = np.array(
            [[0] + list(run.yv), [1] + list(1.0 / (np.exp(run.yv) + 1))]
        )
        run.rho[0, 1, 1] = np.array(
            [[0] + list(run.yv), [1] + list(2.0 / (np.exp(run.yv) + 1))]
        )
        del run.full
        with patch("matplotlib.pyplot.plot") as _plt:
            run.plotRhoOffDiagY(0, 1, 2.5)
            run.plotRhoOffDiagY(0, 1, 2.5, mass=True)
            self.assertEqual(_plt.call_count, 0)
        run.full = True
        with patch("matplotlib.pyplot.plot") as _plt, patch(
            "matplotlib.pyplot.xscale"
        ) as _xs, patch("matplotlib.pyplot.xlabel") as _xl, patch(
            "matplotlib.pyplot.ylabel"
        ) as _yl:
            run.plotRhoOffDiagY(0, 1, 2.5, im=False)
            _xs.assert_called_once_with("log")
            _xl.assert_called_once_with("$x$")
            _yl.assert_called_once_with(r"$\rho_{\alpha\beta}$")
            _plt.assert_called_once()
            xv, yv = run.interpolateRhoIJ(0, 1, 2.5, ri=0, mass=False)
            self.assertEqualArray(_plt.call_args[0], [xv, yv])
            self.assertEqual(
                _plt.call_args[1],
                {
                    "label": r"%s $\alpha\beta$=%d%d re" % (run.label, 1, 2),
                    "ls": "-",
                    "c": "k",
                },
            )
            run.plotRhoOffDiagY(0, 1, 2.5, lc="r")
            self.assertEqual(_plt.call_count, 3)
            xv, yv = run.interpolateRhoIJ(0, 1, 2.5, ri=1, mass=False)
            self.assertEqualArray(_plt.call_args[0], [xv, yv])
            self.assertEqual(
                _plt.call_args[1],
                {
                    "label": r"%s $\alpha\beta$=%d%d im" % (run.label, 1, 2),
                    "ls": ":",
                    "c": "r",
                },
            )
        with patch("matplotlib.pyplot.plot") as _plt:
            run.plotRhoOffDiagY(0, 1, 2.5, lc="r", im=False)
            xv, yv = run.interpolateRhoIJ(0, 1, 2.5, ri=0, mass=False)
            self.assertEqualArray(_plt.call_args[0], [xv, yv])
            self.assertEqual(
                _plt.call_args[1],
                {
                    "label": r"%s $\alpha\beta$=%d%d re" % (run.label, 1, 2),
                    "ls": "-",
                    "c": "r",
                },
            )
        run.rhoM[1, 0, 0] = np.array(
            [
                [0] + list(1.0 / (np.exp(run.yv) + 1)),
                [1] + list(2.0 / (np.exp(run.yv) + 1)),
            ]
        )
        with patch("matplotlib.pyplot.plot") as _plt:
            run.plotRhoOffDiagY(1, 0, 2.5, im=False, mass=True)
            xv, yv = run.interpolateRhoIJ(1, 0, 2.5, ri=0, mass=True)
            self.assertEqualArray(_plt.call_args[0], [xv, yv])

    def test_plotdRhoOffDiagY(self):
        """test plotdRhoOffDiagY"""
        run = fpom.FortEPiaNORun("output/nonexistent")
        run.full = True
        with patch("matplotlib.pyplot.plot") as _plt:
            run.plotdRhoOffDiagY(0, 1, 2.5)
            run.plotdRhoOffDiagY(0, 1, 2.5, mass=True)
            self.assertEqual(_plt.call_count, 0)
        run.rho = np.asarray([])
        run.rhoM = np.asarray([])
        run.rho = np.asarray(
            [
                [
                    [None, None],
                    [None, None],
                ],
                [
                    [None, None],
                    [None, None],
                ],
            ],
        )
        run.rhoM = np.asarray(
            [
                [
                    [None, None],
                    [None, None],
                ],
                [
                    [None, None],
                    [None, None],
                ],
            ],
        )
        run.yv = np.linspace(0.01, 20, 200)
        run.rho[0, 1, 0] = np.array(
            [[0] + list(run.yv), [1] + list(1.0 / (np.exp(run.yv) + 1))]
        )
        run.rho[0, 1, 1] = np.array(
            [[0] + list(run.yv), [1] + list(2.0 / (np.exp(run.yv) + 1))]
        )
        del run.full
        with patch("matplotlib.pyplot.plot") as _plt:
            run.plotdRhoOffDiagY(0, 1, 2.5)
            run.plotdRhoOffDiagY(0, 1, 2.5, mass=True)
            self.assertEqual(_plt.call_count, 0)
        run.full = True
        with patch("matplotlib.pyplot.plot") as _plt, patch(
            "matplotlib.pyplot.xscale"
        ) as _xs, patch("matplotlib.pyplot.xlabel") as _xl, patch(
            "matplotlib.pyplot.ylabel"
        ) as _yl:
            run.plotdRhoOffDiagY(0, 1, 2.5, im=False)
            _xs.assert_called_once_with("log")
            _xl.assert_called_once_with("$x$")
            _yl.assert_called_once_with(r"$d\rho_{\alpha\beta}/dx$")
            _plt.assert_called_once()
            xv, yv = run.interpolateRhoIJ(0, 1, 2.5, ri=0, mass=False)
            self.assertEqualArray(_plt.call_args[0], [xv, np.gradient(yv, xv)])
            self.assertEqual(
                _plt.call_args[1],
                {
                    "label": r"%s $\alpha\beta$=%d%d re" % (run.label, 1, 2),
                    "ls": "-",
                    "c": "k",
                },
            )
            run.plotdRhoOffDiagY(0, 1, 2.5, lc="r")
            self.assertEqual(_plt.call_count, 3)
            xv, yv = run.interpolateRhoIJ(0, 1, 2.5, ri=1, mass=False)
            self.assertEqualArray(_plt.call_args[0], [xv, np.gradient(yv, xv)])
            self.assertEqual(
                _plt.call_args[1],
                {
                    "label": r"%s $\alpha\beta$=%d%d im" % (run.label, 1, 2),
                    "ls": ":",
                    "c": "r",
                },
            )
        with patch("matplotlib.pyplot.plot") as _plt:
            run.plotdRhoOffDiagY(0, 1, 2.5, lc="r", im=False)
            xv, yv = run.interpolateRhoIJ(0, 1, 2.5, ri=0, mass=False)
            self.assertEqualArray(_plt.call_args[0], [xv, np.gradient(yv, xv)])
            self.assertEqual(
                _plt.call_args[1],
                {
                    "label": r"%s $\alpha\beta$=%d%d re" % (run.label, 1, 2),
                    "ls": "-",
                    "c": "r",
                },
            )
        run.rhoM[1, 0, 0] = np.array(
            [
                [0] + list(1.0 / (np.exp(run.yv) + 1)),
                [1] + list(2.0 / (np.exp(run.yv) + 1)),
            ]
        )
        with patch("matplotlib.pyplot.plot") as _plt:
            run.plotdRhoOffDiagY(1, 0, 2.5, im=False, mass=True)
            xv, yv = run.interpolateRhoIJ(1, 0, 2.5, ri=0, mass=True)
            self.assertEqualArray(_plt.call_args[0], [xv, np.gradient(yv, xv)])

    def test_plotNeff(self):
        """test plotNeff"""
        run = fpom.FortEPiaNORun("output/nonexistent")
        run.nnu = 2
        run.zCol = 1
        run.lowReheating = False
        with patch("matplotlib.pyplot.plot") as _plt:
            run.plotNeff()
            run.Neffdat = np.array([[np.nan, np.nan], [np.nan, np.nan]])
            run.plotNeff()
            run.rho = self.explanatory.rho
            run.plotNeff()
            self.assertEqual(_plt.call_count, 0)
        run.zdat = np.array([[0.1, 1.1, 1.0], [1, 1.2, 1.1], [10, 1.3, 1.1]])
        data = [
            [
                0.1,
                8.0 / 7.0 * 3.0 / (fpom.PISQD15 * 1.1 ** 4),
                8.0
                / 7.0
                * 3.0
                / (fpom.PISQD15 * 1.1 ** 4)
                * (11.0 / 4.0) ** (4.0 / 3.0),
            ],
            [
                1.0,
                8.0 / 7.0 * 3.3 / (fpom.PISQD15 * 1.2 ** 4),
                8.0
                / 7.0
                * 3.3
                / (fpom.PISQD15 * 1.2 ** 4)
                * (11.0 / 4.0) ** (4.0 / 3.0),
            ],
            [
                10.0,
                8.0 / 7.0 * 3.6 / (fpom.PISQD15 * 1.3 ** 4),
                8.0
                / 7.0
                * 3.6
                / (fpom.PISQD15 * 1.3 ** 4)
                * (11.0 / 4.0) ** (4.0 / 3.0),
            ],
        ]
        data = np.asarray(data)
        with patch("matplotlib.pyplot.plot") as _plt, patch(
            "numpy.savetxt"
        ) as _sv, patch(
            "fortepianoOutput.FortEPiaNORun.integrateRho_yn",
            side_effect=[1.0, 2.0, 1.1, 2.2, 1.2, 2.4],
        ) as _int:
            run.plotNeff()
            _plt.assert_called_once()
            self.assertEqualArray(
                _plt.call_args[0],
                fpom.stripRepeated(data, 0, 1),
            )
            self.assertEqual(
                _plt.call_args[1],
                {"ls": "-", "c": "k", "label": run.label},
            )
            self.assertEqual(_sv.call_args[0][0], os.path.join(run.folder, "Neff.dat"))
            self.assertEqualArray(_sv.call_args[0][1], data)
            self.assertEqual(_sv.call_args[1], {"fmt": "%.7e"})
            _int.assert_has_calls(
                [call(ii, 3, ix=jj) for ii in [0, 1] for jj in [0, 1, 2]],
                any_order=True,
            )
        del run.rho
        del run.zdat
        run.Neffdat = data
        with patch("matplotlib.pyplot.plot") as _plt, patch(
            "numpy.savetxt"
        ) as _sv, patch("matplotlib.pyplot.xlabel") as _xl, patch(
            "matplotlib.pyplot.xscale"
        ) as _xs, patch(
            "fortepianoOutput.FortEPiaNORun.integrateRho_yn"
        ) as _int:
            run.plotNeff(axes=False, lc="r", ls=":", lab="mylabel")
            _xs.assert_called_once_with("log")
            _xl.assert_called_once_with("$x$")
            _plt.assert_called_once()
            self.assertEqualArray(
                _plt.call_args[0],
                fpom.stripRepeated(data, 0, 1),
            )
            self.assertEqual(
                _plt.call_args[1],
                {"ls": ":", "c": "r", "label": "mylabel"},
            )
            self.assertEqual(_sv.call_count, 0)
            self.assertEqual(_int.call_count, 0)
        plt.figure()
        ax = plt.gca()
        ax1 = ax.twiny()
        with patch("matplotlib.pyplot.plot") as _plt, patch(
            "numpy.savetxt"
        ) as _sv, patch("matplotlib.pyplot.xlabel") as _xl, patch(
            "matplotlib.pyplot.xscale"
        ) as _xs, patch(
            "fortepianoOutput.FortEPiaNORun.integrateRho_yn"
        ) as _int, patch(
            "matplotlib.pyplot.gca", return_value=ax
        ) as _gca, patch(
            "matplotlib.pyplot.Axes.twinx", return_value=ax1
        ) as _twx, patch(
            "matplotlib.pyplot.Axes.set_ylim"
        ) as _syl, patch(
            "matplotlib.pyplot.Axes.set_ylabel"
        ) as _sya:
            run.plotNeff()
            _xs.assert_called_once_with("log")
            _xl.assert_called_once_with("$x$")
            _plt.assert_called_once()
            self.assertEqualArray(
                _plt.call_args[0],
                fpom.stripRepeated(data, 0, 1),
            )
            self.assertEqual(_sv.call_count, 0)
            self.assertEqual(_int.call_count, 0)
            _gca.assert_called_once_with()
            _twx.assert_called_once_with()
            _syl.assert_any_call([0.5, 4.5])
            self.assertEqualArray(
                _syl.call_args_list[1][0][0],
                np.asarray([0.5, 4.5]) * (11.0 / 4) ** (4.0 / 3),
            )
            _sya.assert_any_call(r"$N_{\rm eff}^{\rm in}$")
            _sya.assert_any_call(r"$N_{\rm eff}^{\rm now}$")
            run.plotNeff(nefflims=[1.5, 2.5])
            self.assertEqualArray(_syl.call_args_list[2][0][0], [1.5, 2.5])
            self.assertEqualArray(
                _syl.call_args_list[3][0][0],
                np.asarray([1.5, 2.5]) * (11.0 / 4) ** (4.0 / 3),
            )

        # test with low reheating
        run = fpom.FortEPiaNORun("output/nonexistent")
        run.nnu = 2
        run.zCol = 2
        run.lowReheating = True
        with patch("matplotlib.pyplot.plot") as _plt:
            run.plotNeff()
            run.Neffdat = np.array([[np.nan, np.nan, np.nan], [np.nan, np.nan, np.nan]])
            run.plotNeff()
            run.rho = self.explanatory.rho
            run.plotNeff()
            self.assertEqual(_plt.call_count, 0)
        run.zdat = np.array(
            [[0.1, 0.002, 1.1, 1.0], [1, 0.003, 1.2, 1.1], [10, 0.004, 1.3, 1.1]]
        )
        data = [
            [
                0.1,
                0.002,
                8.0 / 7.0 * 3.0 / (fpom.PISQD15 * 1.1 ** 4),
                8.0
                / 7.0
                * 3.0
                / (fpom.PISQD15 * 1.1 ** 4)
                * (11.0 / 4.0) ** (4.0 / 3.0),
            ],
            [
                1.0,
                0.003,
                8.0 / 7.0 * 3.3 / (fpom.PISQD15 * 1.2 ** 4),
                8.0
                / 7.0
                * 3.3
                / (fpom.PISQD15 * 1.2 ** 4)
                * (11.0 / 4.0) ** (4.0 / 3.0),
            ],
            [
                10.0,
                0.004,
                8.0 / 7.0 * 3.6 / (fpom.PISQD15 * 1.3 ** 4),
                8.0
                / 7.0
                * 3.6
                / (fpom.PISQD15 * 1.3 ** 4)
                * (11.0 / 4.0) ** (4.0 / 3.0),
            ],
        ]
        data = np.asarray(data)
        # test plotting x
        with patch("matplotlib.pyplot.plot") as _plt, patch(
            "numpy.savetxt"
        ) as _sv, patch(
            "fortepianoOutput.FortEPiaNORun.integrateRho_yn",
            side_effect=[1.0, 2.0, 1.1, 2.2, 1.2, 2.4],
        ) as _int:
            run.plotNeff()
            _plt.assert_called_once()
            self.assertEqualArray(
                _plt.call_args[0],
                fpom.stripRepeated(data, 0, 2),
            )
            self.assertEqual(
                _plt.call_args[1],
                {"ls": "-", "c": "k", "label": run.label},
            )
            self.assertEqual(_sv.call_args[0][0], os.path.join(run.folder, "Neff.dat"))
            self.assertEqualArray(_sv.call_args[0][1], data)
            self.assertEqual(_sv.call_args[1], {"fmt": "%.7e"})
            _int.assert_has_calls(
                [call(ii, 3, ix=jj) for ii in [0, 1] for jj in [0, 1, 2]],
                any_order=True,
            )
        # test plotting t
        with patch("matplotlib.pyplot.plot") as _plt, patch(
            "numpy.savetxt"
        ) as _sv, patch(
            "fortepianoOutput.FortEPiaNORun.integrateRho_yn",
            side_effect=[1.0, 2.0, 1.1, 2.2, 1.2, 2.4],
        ) as _int:
            run.plotNeff(useT=True)
            _plt.assert_called_once()
            self.assertEqualArray(
                _plt.call_args[0],
                fpom.stripRepeated(data, 1, 2),
            )
            self.assertEqual(
                _plt.call_args[1],
                {"ls": "-", "c": "k", "label": run.label},
            )
            self.assertEqual(_sv.call_args[0][0], os.path.join(run.folder, "Neff.dat"))
            self.assertEqualArray(_sv.call_args[0][1], data)
            self.assertEqual(_sv.call_args[1], {"fmt": "%.7e"})
            _int.assert_has_calls(
                [call(ii, 3, ix=jj) for ii in [0, 1] for jj in [0, 1, 2]],
                any_order=True,
            )
        del run.rho
        del run.zdat
        run.Neffdat = data
        with patch("matplotlib.pyplot.plot") as _plt, patch(
            "numpy.savetxt"
        ) as _sv, patch("matplotlib.pyplot.xlabel") as _xl, patch(
            "matplotlib.pyplot.xscale"
        ) as _xs, patch(
            "fortepianoOutput.FortEPiaNORun.integrateRho_yn"
        ) as _int:
            run.plotNeff(axes=False, lc="r", ls=":", lab="mylabel")
            _xs.assert_called_once_with("log")
            _xl.assert_called_once_with("$x$")
            _plt.assert_called_once()
            self.assertEqualArray(
                _plt.call_args[0],
                fpom.stripRepeated(data, 0, 2),
            )
            self.assertEqual(
                _plt.call_args[1],
                {"ls": ":", "c": "r", "label": "mylabel"},
            )
            self.assertEqual(_sv.call_count, 0)
            self.assertEqual(_int.call_count, 0)
        plt.figure()
        ax = plt.gca()
        ax1 = ax.twiny()
        with patch("matplotlib.pyplot.plot") as _plt, patch(
            "numpy.savetxt"
        ) as _sv, patch("matplotlib.pyplot.xlabel") as _xl, patch(
            "matplotlib.pyplot.xscale"
        ) as _xs, patch(
            "fortepianoOutput.FortEPiaNORun.integrateRho_yn"
        ) as _int, patch(
            "matplotlib.pyplot.gca", return_value=ax
        ) as _gca, patch(
            "matplotlib.pyplot.Axes.twinx", return_value=ax1
        ) as _twx, patch(
            "matplotlib.pyplot.Axes.set_ylim"
        ) as _syl, patch(
            "matplotlib.pyplot.Axes.set_ylabel"
        ) as _sya:
            run.plotNeff()
            _xs.assert_called_once_with("log")
            _xl.assert_called_once_with("$x$")
            _plt.assert_called_once()
            self.assertEqualArray(
                _plt.call_args[0],
                fpom.stripRepeated(data, 0, 2),
            )
            self.assertEqual(_sv.call_count, 0)
            self.assertEqual(_int.call_count, 0)
            _gca.assert_called_once_with()
            _twx.assert_called_once_with()
            _syl.assert_any_call([0.5, 4.5])
            self.assertEqualArray(
                _syl.call_args_list[1][0][0],
                np.asarray([0.5, 4.5]) * (11.0 / 4) ** (4.0 / 3),
            )
            _sya.assert_any_call(r"$N_{\rm eff}^{\rm in}$")
            _sya.assert_any_call(r"$N_{\rm eff}^{\rm now}$")
            run.plotNeff(nefflims=[1.5, 2.5])
            self.assertEqualArray(_syl.call_args_list[2][0][0], [1.5, 2.5])
            self.assertEqualArray(
                _syl.call_args_list[3][0][0],
                np.asarray([1.5, 2.5]) * (11.0 / 4) ** (4.0 / 3),
            )

    def test_plotEnergyDensity(self):
        """test plotEnergyDensity"""
        run = fpom.FortEPiaNORun("output/nonexistent")
        with patch("matplotlib.pyplot.plot") as _plt:
            run.plotEnergyDensity()
            run.endens = None
            run.plotEnergyDensity()
            run.endens = [0, 1, 2]
            run.plotEnergyDensity()
            self.assertEqual(_plt.call_count, 0)
        run = self.explanatory
        colors = ["r", "b", "g", "#ff9933", "#ff9933", "#ff9933", "#ff00ff"]
        styles = ["-", "-", "-", ":", "-.", "--", "-"]
        skip = [False, False, False, False, False, False, False]
        with patch("matplotlib.pyplot.plot") as _plt:
            run.plotEnergyDensity()
            self.assertEqual(_plt.call_count, 10 if run.lowReheating else 9)
            self.assertEqualArray(
                _plt.call_args_list[0][0],
                [
                    run.endens[:, 0],
                    np.asarray([np.sum(cl[run.zCol + 1 :]) for cl in run.endens]),
                ],
            )
            self.assertEqual(
                _plt.call_args_list[0][1],
                {
                    "label": "total",
                    "c": "k",
                    "ls": "-",
                    "lw": 1,
                },
            )
            for ix, lab in enumerate(
                [
                    r"$\gamma$",
                    "$e$",
                    r"$\mu$",
                    r"$\nu_e$",
                    r"$\nu_\mu$",
                    r"$\nu_\tau$",
                ]
            ):
                self.assertEqualArray(
                    _plt.call_args_list[1 + ix][0],
                    [run.endens[:, 0], run.endens[:, run.zCol + 1 + ix]],
                )
                self.assertEqual(
                    _plt.call_args_list[1 + ix][1],
                    {
                        "label": lab,
                        "c": colors[ix],
                        "ls": styles[ix],
                        "lw": 1,
                    },
                )
            if run.lowReheating:
                self.assertEqualArray(
                    _plt.call_args_list[2 + ix][0],
                    [run.endens[:, 0], run.endens[:, run.zCol + 2 + ix]],
                )
                self.assertEqual(
                    _plt.call_args_list[2 + ix][1],
                    {
                        "label": r"$\phi$",
                        "c": colors[ix + 1],
                        "ls": styles[ix + 1],
                        "lw": 1,
                    },
                )
            self.assertEqualArray(
                _plt.call_args_list[-2][0],
                [
                    run.endens[:, 0],
                    run.endens[:, run.zCol + 1] + run.endens[:, run.zCol + 2],
                ],
            )
            self.assertEqual(
                _plt.call_args_list[-2][1],
                {
                    "label": r"$\gamma+e$",
                    "c": "#00ccff",
                    "ls": ":",
                    "lw": 1,
                },
            )
            self.assertEqualArray(
                _plt.call_args_list[-1][0],
                [
                    run.endens[:, 0],
                    run.endens[:, run.zCol + 1]
                    + run.endens[:, run.zCol + 2]
                    + run.endens[:, run.zCol + 3],
                ],
            )
            self.assertEqual(
                _plt.call_args_list[-1][1],
                {
                    "label": r"$\gamma+e+\mu$",
                    "c": "#6666ff",
                    "ls": "--",
                    "lw": 1,
                },
            )
        # some tweaks on the inputs
        run.endens = np.c_[run.endens, run.endens[:, 7]]
        colors = ["b", "w", "m", "#ff9933", "c", "#ff9933", "#ff00ff"]
        styles = [":", "-", ".", ":", ".", ":", "-"]
        skip = [False, True, False, False, True, True, False]
        labels = [
            r"a$\gamma$",
            "a$e$",
            r"a$\mu$",
            r"a$\nu_e$",
            r"a$\nu_\mu$",
            r"a$\nu_\tau$",
            r"a$\nu_s$",
        ]
        with patch("matplotlib.pyplot.plot") as _plt:
            run.plotEnergyDensity(
                gamma_e=True,
                gec="#002233",
                ges="-",
                gamma_e_mu=True,
                gemc="#557788",
                gems=":",
                labels=labels,
                colors=colors,
                styles=styles,
                skip=skip,
                lw=2,
            )
            self.assertEqual(_plt.call_count, 7)
            self.assertEqualArray(
                _plt.call_args_list[0][0],
                [
                    run.endens[:, 0],
                    np.asarray([np.sum(cl[run.zCol + 1 :]) for cl in run.endens]),
                ],
            )
            self.assertEqual(
                _plt.call_args_list[0][1],
                {
                    "label": "total",
                    "c": "k",
                    "ls": "-",
                    "lw": 2,
                },
            )
            ii = 1
            for ix, lab in enumerate(labels):
                if skip[ix]:
                    continue
                self.assertEqualArray(
                    _plt.call_args_list[ii][0],
                    [run.endens[:, 0], run.endens[:, run.zCol + 1 + ix]],
                )
                self.assertEqual(
                    _plt.call_args_list[ii][1],
                    {
                        "label": lab,
                        "c": colors[ix],
                        "ls": styles[ix],
                        "lw": 2,
                    },
                )
                ii += 1
            self.assertEqualArray(
                _plt.call_args_list[-2][0],
                [
                    run.endens[:, 0],
                    run.endens[:, run.zCol + 1] + run.endens[:, run.zCol + 2],
                ],
            )
            self.assertEqual(
                _plt.call_args_list[-2][1],
                {
                    "label": r"$\gamma+e$",
                    "c": "#002233",
                    "ls": "-",
                    "lw": 2,
                },
            )
            self.assertEqualArray(
                _plt.call_args_list[-1][0],
                [
                    run.endens[:, 0],
                    run.endens[:, run.zCol + 1]
                    + run.endens[:, run.zCol + 2]
                    + run.endens[:, run.zCol + 3],
                ],
            )
            self.assertEqual(
                _plt.call_args_list[-1][1],
                {
                    "label": r"$\gamma+e+\mu$",
                    "c": "#557788",
                    "ls": ":",
                    "lw": 2,
                },
            )
        skip = [False, False, False, False, False, False, False]
        with patch("matplotlib.pyplot.plot") as _plt:
            run.plotEnergyDensity(
                gamma_e=False,
                gec="#002233",
                ges="-",
                gamma_e_mu=False,
                gemc="#557788",
                gems=":",
                labels=labels,
                colors=colors,
                styles=styles,
                skip=skip,
                lw=2,
            )
            self.assertEqual(_plt.call_count, 8)
            self.assertEqualArray(
                _plt.call_args_list[0][0],
                [
                    run.endens[:, 0],
                    np.asarray([np.sum(cl[run.zCol + 1 :]) for cl in run.endens]),
                ],
            )
            self.assertEqual(
                _plt.call_args_list[0][1],
                {
                    "label": "total",
                    "c": "k",
                    "ls": "-",
                    "lw": 2,
                },
            )
            ii = 1
            for ix, lab in enumerate(labels):
                if skip[ix]:
                    continue
                self.assertEqualArray(
                    _plt.call_args_list[ii][0],
                    [run.endens[:, 0], run.endens[:, run.zCol + 1 + ix]],
                )
                self.assertEqual(
                    _plt.call_args_list[ii][1],
                    {
                        "label": lab,
                        "c": colors[ix],
                        "ls": styles[ix],
                        "lw": 2,
                    },
                )
                ii += 1
        with patch("matplotlib.pyplot.plot") as _plt:
            run.plotEnergyDensity(
                gamma_e=True,
                gec="#002233",
                ges="-",
                gamma_e_mu=True,
                gemc="#557788",
                gems=":",
                labels=labels,
                colors=colors,
                styles=styles,
                skip=skip,
                allstyles="--",
                alllabels="newlab",
                lw=2,
            )
            self.assertEqual(_plt.call_count, 10)
            self.assertEqualArray(
                _plt.call_args_list[0][0],
                [
                    run.endens[:, 0],
                    np.asarray([np.sum(cl[run.zCol + 1 :]) for cl in run.endens]),
                ],
            )
            self.assertEqual(
                _plt.call_args_list[0][1],
                {
                    "label": "newlab",
                    "c": "k",
                    "ls": "--",
                    "lw": 2,
                },
            )
            ii = 1
            for ix, lab in enumerate(labels):
                if skip[ix]:
                    continue
                self.assertEqualArray(
                    _plt.call_args_list[ii][0],
                    [run.endens[:, 0], run.endens[:, run.zCol + 1 + ix]],
                )
                self.assertEqual(
                    _plt.call_args_list[ii][1],
                    {
                        "label": "newlab",
                        "c": colors[ix],
                        "ls": "--",
                        "lw": 2,
                    },
                )
                ii += 1
            self.assertEqualArray(
                _plt.call_args_list[-2][0],
                [
                    run.endens[:, 0],
                    run.endens[:, run.zCol + 1] + run.endens[:, run.zCol + 2],
                ],
            )
            self.assertEqual(
                _plt.call_args_list[-2][1],
                {
                    "label": "newlab",
                    "c": "#002233",
                    "ls": "--",
                    "lw": 2,
                },
            )
            self.assertEqualArray(
                _plt.call_args_list[-1][0],
                [
                    run.endens[:, 0],
                    run.endens[:, run.zCol + 1]
                    + run.endens[:, run.zCol + 2]
                    + run.endens[:, run.zCol + 3],
                ],
            )
            self.assertEqual(
                _plt.call_args_list[-1][1],
                {
                    "label": "newlab",
                    "c": "#557788",
                    "ls": "--",
                    "lw": 2,
                },
            )

    def test_plotEntropy(self):
        """test plotEntropy"""
        run = fpom.FortEPiaNORun("output/nonexistent")
        with patch("matplotlib.pyplot.plot") as _plt:
            run.plotEntropy()
            run.entropy = None
            run.plotEntropy()
            run.entropy = [0, 1, 2]
            run.plotEntropy()
            self.assertEqual(_plt.call_count, 0)
        run = self.explanatory
        colors = ["r", "b", "g", "#ff9933", "#ff9933", "#ff9933", "#ff00ff"]
        styles = ["-", "-", "-", ":", "-.", "--", "-"]
        skip = [False, False, False, False, False, False, False]
        with patch("matplotlib.pyplot.plot") as _plt:
            run.plotEntropy()
            self.assertEqual(_plt.call_count, 9)
            self.assertEqualArray(
                _plt.call_args_list[0][0],
                [
                    run.entropy[:, 0],
                    np.asarray([np.sum(cl[run.zCol + 1 :]) for cl in run.entropy]),
                ],
            )
            self.assertEqual(
                _plt.call_args_list[0][1],
                {
                    "label": "total",
                    "c": "k",
                    "ls": "-",
                    "lw": 1,
                },
            )
            for ix, lab in enumerate(
                [
                    r"$\gamma$",
                    "$e$",
                    r"$\mu$",
                    r"$\nu_e$",
                    r"$\nu_\mu$",
                    r"$\nu_\tau$",
                ]
            ):
                self.assertEqualArray(
                    _plt.call_args_list[1 + ix][0],
                    [run.entropy[:, 0], run.entropy[:, run.zCol + 1 + ix]],
                )
                self.assertEqual(
                    _plt.call_args_list[1 + ix][1],
                    {
                        "label": lab,
                        "c": colors[ix],
                        "ls": styles[ix],
                        "lw": 1,
                    },
                )
            self.assertEqualArray(
                _plt.call_args_list[-2][0],
                [
                    run.entropy[:, 0],
                    run.entropy[:, run.zCol + 1] + run.entropy[:, run.zCol + 2],
                ],
            )
            self.assertEqual(
                _plt.call_args_list[-2][1],
                {
                    "label": r"$\gamma+e$",
                    "c": "#00ccff",
                    "ls": ":",
                    "lw": 1,
                },
            )
            self.assertEqualArray(
                _plt.call_args_list[-1][0],
                [
                    run.entropy[:, 0],
                    run.entropy[:, run.zCol + 1]
                    + run.entropy[:, run.zCol + 2]
                    + run.entropy[:, run.zCol + 3],
                ],
            )
            self.assertEqual(
                _plt.call_args_list[-1][1],
                {
                    "label": r"$\gamma+e+\mu$",
                    "c": "#6666ff",
                    "ls": "--",
                    "lw": 1,
                },
            )
        # some tweaks on the inputs
        run.entropy = np.c_[run.entropy, run.entropy[:, 7]]
        colors = ["b", "w", "m", "#ff9933", "c", "#ff9933", "#ff00ff"]
        styles = [":", "-", ".", ":", ".", ":", "-"]
        skip = [False, True, False, False, True, True, False]
        labels = [
            r"a$\gamma$",
            "a$e$",
            r"a$\mu$",
            r"a$\nu_e$",
            r"a$\nu_\mu$",
            r"a$\nu_\tau$",
            r"a$\nu_s$",
        ]
        with patch("matplotlib.pyplot.plot") as _plt:
            run.plotEntropy(
                gamma_e=True,
                gec="#002233",
                ges="-",
                gamma_e_mu=True,
                gemc="#557788",
                gems=":",
                labels=labels,
                colors=colors,
                styles=styles,
                skip=skip,
                lw=2,
            )
            self.assertEqual(_plt.call_count, 7)
            self.assertEqualArray(
                _plt.call_args_list[0][0],
                [
                    run.entropy[:, 0],
                    np.asarray([np.sum(cl[run.zCol + 1 :]) for cl in run.entropy]),
                ],
            )
            self.assertEqual(
                _plt.call_args_list[0][1],
                {
                    "label": "total",
                    "c": "k",
                    "ls": "-",
                    "lw": 2,
                },
            )
            ii = 1
            for ix, lab in enumerate(labels):
                if skip[ix]:
                    continue
                self.assertEqualArray(
                    _plt.call_args_list[ii][0],
                    [run.entropy[:, 0], run.entropy[:, run.zCol + 1 + ix]],
                )
                self.assertEqual(
                    _plt.call_args_list[ii][1],
                    {
                        "label": lab,
                        "c": colors[ix],
                        "ls": styles[ix],
                        "lw": 2,
                    },
                )
                ii += 1
            self.assertEqualArray(
                _plt.call_args_list[-2][0],
                [
                    run.entropy[:, 0],
                    run.entropy[:, run.zCol + 1] + run.entropy[:, run.zCol + 2],
                ],
            )
            self.assertEqual(
                _plt.call_args_list[-2][1],
                {
                    "label": r"$\gamma+e$",
                    "c": "#002233",
                    "ls": "-",
                    "lw": 2,
                },
            )
            self.assertEqualArray(
                _plt.call_args_list[-1][0],
                [
                    run.entropy[:, 0],
                    run.entropy[:, run.zCol + 1]
                    + run.entropy[:, run.zCol + 2]
                    + run.entropy[:, run.zCol + 3],
                ],
            )
            self.assertEqual(
                _plt.call_args_list[-1][1],
                {
                    "label": r"$\gamma+e+\mu$",
                    "c": "#557788",
                    "ls": ":",
                    "lw": 2,
                },
            )
        skip = [False, False, False, False, False, False, False]
        with patch("matplotlib.pyplot.plot") as _plt:
            run.plotEntropy(
                gamma_e=False,
                gec="#002233",
                ges="-",
                gamma_e_mu=False,
                gemc="#557788",
                gems=":",
                labels=labels,
                colors=colors,
                styles=styles,
                skip=skip,
                lw=2,
            )
            self.assertEqual(_plt.call_count, 8)
            self.assertEqualArray(
                _plt.call_args_list[0][0],
                [
                    run.entropy[:, 0],
                    np.asarray([np.sum(cl[run.zCol + 1 :]) for cl in run.entropy]),
                ],
            )
            self.assertEqual(
                _plt.call_args_list[0][1],
                {
                    "label": "total",
                    "c": "k",
                    "ls": "-",
                    "lw": 2,
                },
            )
            ii = 1
            for ix, lab in enumerate(labels):
                if skip[ix]:
                    continue
                self.assertEqualArray(
                    _plt.call_args_list[ii][0],
                    [run.entropy[:, 0], run.entropy[:, run.zCol + 1 + ix]],
                )
                self.assertEqual(
                    _plt.call_args_list[ii][1],
                    {
                        "label": lab,
                        "c": colors[ix],
                        "ls": styles[ix],
                        "lw": 2,
                    },
                )
                ii += 1
        with patch("matplotlib.pyplot.plot") as _plt:
            run.plotEntropy(
                gamma_e=True,
                gec="#002233",
                ges="-",
                gamma_e_mu=True,
                gemc="#557788",
                gems=":",
                labels=labels,
                colors=colors,
                styles=styles,
                skip=skip,
                allstyles="--",
                alllabels="newlab",
                lw=2,
            )
            self.assertEqual(_plt.call_count, 10)
            self.assertEqualArray(
                _plt.call_args_list[0][0],
                [
                    run.entropy[:, 0],
                    np.asarray([np.sum(cl[run.zCol + 1 :]) for cl in run.entropy]),
                ],
            )
            self.assertEqual(
                _plt.call_args_list[0][1],
                {
                    "label": "newlab",
                    "c": "k",
                    "ls": "--",
                    "lw": 2,
                },
            )
            ii = 1
            for ix, lab in enumerate(labels):
                if skip[ix]:
                    continue
                self.assertEqualArray(
                    _plt.call_args_list[ii][0],
                    [run.entropy[:, 0], run.entropy[:, run.zCol + 1 + ix]],
                )
                self.assertEqual(
                    _plt.call_args_list[ii][1],
                    {
                        "label": "newlab",
                        "c": colors[ix],
                        "ls": "--",
                        "lw": 2,
                    },
                )
                ii += 1
            self.assertEqualArray(
                _plt.call_args_list[-2][0],
                [
                    run.entropy[:, 0],
                    run.entropy[:, run.zCol + 1] + run.entropy[:, run.zCol + 2],
                ],
            )
            self.assertEqual(
                _plt.call_args_list[-2][1],
                {
                    "label": "newlab",
                    "c": "#002233",
                    "ls": "--",
                    "lw": 2,
                },
            )
            self.assertEqualArray(
                _plt.call_args_list[-1][0],
                [
                    run.entropy[:, 0],
                    run.entropy[:, run.zCol + 1]
                    + run.entropy[:, run.zCol + 2]
                    + run.entropy[:, run.zCol + 3],
                ],
            )
            self.assertEqual(
                _plt.call_args_list[-1][1],
                {
                    "label": "newlab",
                    "c": "#557788",
                    "ls": "--",
                    "lw": 2,
                },
            )

    def test_plotNumberDensity(self):
        """test plotNumberDensity"""
        run = fpom.FortEPiaNORun("output/nonexistent")
        with patch("matplotlib.pyplot.plot") as _plt:
            run.plotNumberDensity()
            run.number = None
            run.plotNumberDensity()
            run.number = [0, 1, 2]
            run.plotNumberDensity()
            self.assertEqual(_plt.call_count, 0)
        run = self.explanatory
        colors = ["r", "b", "g", "#ff9933", "#ff9933", "#ff9933", "#ff00ff"]
        styles = ["-", "-", "-", ":", "-.", "--", "-"]
        skip = [False, False, False, False, False, False, False]
        with patch("matplotlib.pyplot.plot") as _plt:
            run.plotNumberDensity()
            self.assertEqual(_plt.call_count, 6)
            for ix, lab in enumerate(
                [
                    r"$\gamma$",
                    "$e$",
                    r"$\mu$",
                    r"$\nu_e$",
                    r"$\nu_\mu$",
                    r"$\nu_\tau$",
                ]
            ):
                self.assertEqualArray(
                    _plt.call_args_list[ix][0],
                    [run.number[:, 0], run.number[:, run.zCol + 1 + ix]],
                )
                self.assertEqual(
                    _plt.call_args_list[ix][1],
                    {
                        "label": lab,
                        "c": colors[ix],
                        "ls": styles[ix],
                        "lw": 1,
                    },
                )
        # some tweaks on the inputs
        run.number = np.c_[run.number, run.number[:, 7]]
        colors = ["b", "w", "m", "#ff9933", "c", "#ff9933", "#ff00ff"]
        styles = [":", "-", ".", ":", ".", ":", "-"]
        skip = [False, True, False, False, True, True, False]
        labels = [
            r"a$\gamma$",
            "a$e$",
            r"a$\mu$",
            r"a$\nu_e$",
            r"a$\nu_\mu$",
            r"a$\nu_\tau$",
            r"a$\nu_s$",
        ]
        with patch("matplotlib.pyplot.plot") as _plt:
            run.plotNumberDensity(
                labels=labels,
                colors=colors,
                styles=styles,
                skip=skip,
                lw=2,
            )
            self.assertEqual(_plt.call_count, 4)
            ii = 0
            for ix, lab in enumerate(labels):
                if skip[ix]:
                    continue
                self.assertEqualArray(
                    _plt.call_args_list[ii][0],
                    [run.number[:, 0], run.number[:, run.zCol + 1 + ix]],
                )
                self.assertEqual(
                    _plt.call_args_list[ii][1],
                    {
                        "label": lab,
                        "c": colors[ix],
                        "ls": styles[ix],
                        "lw": 2,
                    },
                )
                ii += 1

    def test_plotPArthENoPE(self):
        """test plotPArthENoPE"""
        run = self.explanatory
        with patch("matplotlib.pyplot.plot") as _p:
            with self.assertRaises(AttributeError):
                run.plotPArthENoPE(x="a")
            with self.assertRaises(AttributeError):
                run.plotPArthENoPE(y="a")
            _p.assert_not_called()
            run.plotPArthENoPE()
            _p.assert_called_once()
            self.assertEqualArray(_p.call_args[0][0], run.x[run.filter99])
            self.assertEqualArray(_p.call_args[0][1], run.drhonu_dx[run.filter99])
            self.assertEqual(_p.call_args[1], {})
        for x in (
            ["x", "z", "Tgamma", "t"] if run.lowReheating else ["x", "z", "Tgamma"]
        ):
            with patch("matplotlib.pyplot.plot") as _p:
                run.plotPArthENoPE(x=x)
                self.assertEqualArray(_p.call_args[0][0], getattr(run, x)[run.filter99])
        for y in ["rhonu", "drhonu_dx", "drhonu_dx_savgol", "N_func", "N_savgol"]:
            with patch("matplotlib.pyplot.plot") as _p:
                run.plotPArthENoPE(y=y)
                self.assertEqualArray(_p.call_args[0][1], getattr(run, y)[run.filter99])
        with patch("matplotlib.pyplot.plot") as _p:
            run.plotPArthENoPE(c="k", ls="abc")
            self.assertEqual(_p.call_args[1], {"c": "k", "ls": "abc"})
        with patch("matplotlib.pyplot.plot") as _p:
            run.plotPArthENoPE(x="z", c="k", y="N_func", ls="abc")
            self.assertEqualArray(_p.call_args[0][0], run.z[run.filter99])
            self.assertEqualArray(_p.call_args[0][1], run.N_func[run.filter99])
            self.assertEqual(_p.call_args[1], {"c": "k", "ls": "abc"})
        with patch("matplotlib.pyplot.plot") as _p:
            run.plotPArthENoPE(x="z", c="k", y="N_func", ls="abc", filter99=False)
            self.assertEqualArray(_p.call_args[0][0], run.z)
            self.assertEqualArray(_p.call_args[0][1], run.N_func)
            self.assertEqual(_p.call_args[1], {"c": "k", "ls": "abc"})

    def test_doAllPlots(self):
        """test doAllPlots"""
        run = self.explanatory
        with patch("matplotlib.pyplot.close") as _cl, patch(
            "fortepianoOutput.finalizePlot"
        ) as _fp, patch("fortepianoOutput.FortEPiaNORun.plotZ") as _plZ, patch(
            "fortepianoOutput.FortEPiaNORun.plotW"
        ) as _plW, patch(
            "fortepianoOutput.FortEPiaNORun.plotNeff"
        ) as _plNf, patch(
            "fortepianoOutput.FortEPiaNORun.plotRhoDiagY"
        ) as _plRdy, patch(
            "fortepianoOutput.FortEPiaNORun.plotdRhoDiagY"
        ) as _plDRdy, patch(
            "fortepianoOutput.FortEPiaNORun.plotRhoFin"
        ) as _plRf, patch(
            "fortepianoOutput.FortEPiaNORun.plotRhoOffDiagY"
        ) as _plRoy, patch(
            "fortepianoOutput.FortEPiaNORun.plotdRhoOffDiagY"
        ) as _plDRoy:
            run.doAllPlots()
            _cl.assert_called_once_with()
            _plZ.assert_called_once_with(lc="k", lab="z")
            _plW.assert_called_once_with(lc="k", ls=":", lab="w")
            _plNf.assert_called_once_with(lc="k", axes=False)
            self.assertEqual(_plRdy.call_count, 6)
            self.assertEqual(_plDRdy.call_count, 6)
            for i in range(run.nnu):
                _plRdy.assert_any_call(i, 5.0, fpom.styles[i], lc=fpom.colors[i])
                _plRdy.assert_any_call(
                    i, 5.0, fpom.styles[i], lc=fpom.colors[i], mass=True
                )
                _plDRdy.assert_any_call(i, 5.0, fpom.styles[i], lc=fpom.colors[i])
                _plDRdy.assert_any_call(
                    i, 5.0, fpom.styles[i], lc=fpom.colors[i], mass=True
                )
            self.assertEqual(_plRf.call_count, 6)
            for i in range(run.nnu):
                _plRf.assert_any_call(i, ls=fpom.styles[i], lc=fpom.colors[i], y2=True)
                _plRf.assert_any_call(
                    i, ls=fpom.styles[i], lc=fpom.colors[i], y2=True, mass=True
                )
            self.assertEqual(_plRoy.call_count, 3)
            self.assertEqual(_plDRoy.call_count, 3)
            for i in range(run.nnu):
                for j in range(i + 1, run.nnu):
                    _plRoy.assert_any_call(i, j, 5.0, lc=fpom.colors[2 * i + j - 1])
                    _plDRoy.assert_any_call(i, j, 5.0, lc=fpom.colors[2 * i + j - 1])
            _fp.assert_has_calls(
                [
                    call(
                        "%s/z.pdf" % run.folder, xlab="$x$", ylab=r"$z$", xscale="log"
                    ),
                    call(
                        "%s/rho_diag.pdf" % run.folder,
                        yscale="log",
                    ),
                    call(
                        "%s/drho_diag.pdf" % run.folder,
                    ),
                    call(
                        "%s/rho_mass_diag.pdf" % run.folder,
                        yscale="log",
                    ),
                    call(
                        "%s/drho_mass_diag.pdf" % run.folder,
                    ),
                    call(
                        "%s/rhofin_diag.pdf" % run.folder, xscale="linear", yscale="log"
                    ),
                    call(
                        "%s/rhofin_mass_diag.pdf" % run.folder,
                        xscale="linear",
                        yscale="log",
                    ),
                    call("%s/rho_offdiag.pdf" % run.folder),
                    call("%s/drho_offdiag.pdf" % run.folder),
                    call("%s/Neff.pdf" % run.folder, legend=False, Neff_axes=True),
                ]
            )

        with patch("matplotlib.pyplot.close") as _cl, patch(
            "fortepianoOutput.finalizePlot"
        ) as _fp, patch("fortepianoOutput.FortEPiaNORun.plotZ") as _plZ, patch(
            "fortepianoOutput.FortEPiaNORun.plotW"
        ) as _plW, patch(
            "fortepianoOutput.FortEPiaNORun.plotNeff"
        ) as _plNf, patch(
            "fortepianoOutput.FortEPiaNORun.plotRhoDiagY"
        ) as _plRdy, patch(
            "fortepianoOutput.FortEPiaNORun.plotdRhoDiagY"
        ) as _plDRdy, patch(
            "fortepianoOutput.FortEPiaNORun.plotRhoFin"
        ) as _plRf, patch(
            "fortepianoOutput.FortEPiaNORun.plotRhoOffDiagY"
        ) as _plRoy, patch(
            "fortepianoOutput.FortEPiaNORun.plotdRhoOffDiagY"
        ) as _plDRoy:
            run.doAllPlots(yref=2.5, color="abc")
            _cl.assert_called_once_with()
            _plZ.assert_called_once_with(lc="abc", lab="z")
            _plW.assert_called_once_with(lc="abc", ls=":", lab="w")
            _plNf.assert_called_once_with(lc="abc", axes=False)
            self.assertEqual(_plRdy.call_count, 6)
            self.assertEqual(_plDRdy.call_count, 6)
            for i in range(run.nnu):
                _plRdy.assert_any_call(i, 2.5, fpom.styles[i], lc=fpom.colors[i])
                _plRdy.assert_any_call(
                    i, 2.5, fpom.styles[i], lc=fpom.colors[i], mass=True
                )
                _plDRdy.assert_any_call(i, 2.5, fpom.styles[i], lc=fpom.colors[i])
                _plDRdy.assert_any_call(
                    i, 2.5, fpom.styles[i], lc=fpom.colors[i], mass=True
                )
            self.assertEqual(_plRf.call_count, 6)
            for i in range(run.nnu):
                _plRf.assert_any_call(i, ls=fpom.styles[i], lc=fpom.colors[i], y2=True)
                _plRf.assert_any_call(
                    i, ls=fpom.styles[i], lc=fpom.colors[i], y2=True, mass=True
                )
            self.assertEqual(_plRoy.call_count, 3)
            self.assertEqual(_plDRoy.call_count, 3)
            for i in range(run.nnu):
                for j in range(i + 1, run.nnu):
                    _plRoy.assert_any_call(i, j, 2.5, lc=fpom.colors[2 * i + j - 1])
                    _plDRoy.assert_any_call(i, j, 2.5, lc=fpom.colors[2 * i + j - 1])
            _fp.assert_has_calls(
                [
                    call(
                        "%s/z.pdf" % run.folder, xlab="$x$", ylab=r"$z$", xscale="log"
                    ),
                    call(
                        "%s/rho_diag.pdf" % run.folder,
                        yscale="log",
                    ),
                    call(
                        "%s/drho_diag.pdf" % run.folder,
                    ),
                    call(
                        "%s/rho_mass_diag.pdf" % run.folder,
                        yscale="log",
                    ),
                    call(
                        "%s/drho_mass_diag.pdf" % run.folder,
                    ),
                    call(
                        "%s/rhofin_diag.pdf" % run.folder, xscale="linear", yscale="log"
                    ),
                    call(
                        "%s/rhofin_mass_diag.pdf" % run.folder,
                        xscale="linear",
                        yscale="log",
                    ),
                    call("%s/rho_offdiag.pdf" % run.folder),
                    call("%s/drho_offdiag.pdf" % run.folder),
                    call("%s/Neff.pdf" % run.folder, legend=False, Neff_axes=True),
                ]
            )

    def test_integrateRho_yn(self):
        """test integrateRho_yn"""
        run = fpom.FortEPiaNORun("output/nonexistent")
        with self.assertRaises(AttributeError):
            run.integrateRho_yn(0, 2, show=True)
        with self.assertRaises(AttributeError):
            run.integrateRho_yn(0, 2, mass=True)
        run.rho = np.asarray([])
        run.rhoM = np.asarray([])
        with self.assertRaises(AttributeError):
            run.integrateRho_yn(0, 2)
        run.yv = np.array([x for x in range(3)])
        with self.assertRaises(IndexError):
            run.integrateRho_yn(0, 2)
        run.rho = np.asarray(
            [
                [
                    [None, None],
                    [None, None],
                ],
                [
                    [None, None],
                    [None, None],
                ],
            ],
        )
        run.rhoM = np.asarray(
            [
                [
                    [None, None],
                    [None, None],
                ],
                [
                    [None, None],
                    [None, None],
                ],
            ],
        )
        run.yv = np.linspace(0.01, 20, 200)
        run.rho[0, 0, 0] = np.array(
            [[0] + list(run.yv), [1] + list(1.0 / (np.exp(run.yv) + 1))]
        )
        run.rho[1, 1, 0] = np.array(
            [[0] + list(run.yv), [1] + list(2.0 / (np.exp(run.yv) + 1))]
        )
        self.assertTrue(np.isclose(run.integrateRho_yn(0, 2), 0.182691))
        self.assertTrue(np.isclose(run.integrateRho_yn(1, 1), 0.166662))
        run.rhoM[0, 0, 0] = np.array(
            [
                [0] + list(1.0 / (np.exp(run.yv) + 1)),
                [1] + list(2.0 / (np.exp(run.yv) + 1)),
            ]
        )
        self.assertTrue(np.isclose(run.integrateRho_yn(0, 1, mass=True), 0.166662))
        self.assertTrue(
            np.isclose(run.integrateRho_yn(0, 2, ix=0, mass=True), 0.182691)
        )


class TestPrepareIni(unittest.TestCase):
    """Testing the prepareIni module """

    testIni = "test_ini_file.ini"

    def tearDown(self):
        """remove test files"""
        if os.path.exists(self.testIni):
            os.remove(self.testIni)

    def test_module(self):
        """test module variables"""
        self.assertEqual(pim.params_active, ["dm21", "th12", "dm31", "th13", "th23"])
        self.assertEqual(pim.params_sterile, ["dm41", "th14", "th24", "th34"])
        self.assertIsInstance(pim.default_osc_act, dict)
        for v in pim.default_osc_act.values():
            self.assertIsInstance(v, dict)
            self.assertIn("no", v.keys())
            self.assertIn("io", v.keys())
            for s in v.values():
                self.assertIsInstance(s, dict)
                for p in pim.params_active:
                    self.assertIn(p, s.keys())
        self.assertIsInstance(pim.default_osc_ster, dict)
        for v in pim.default_osc_ster.values():
            self.assertIsInstance(v, dict)
            for p in pim.params_sterile:
                self.assertIn(p, v.keys())

    def test_setParser(self):
        """test setParser"""
        parser = argparse.ArgumentParser()
        parser.add_argument = MagicMock()
        with patch("argparse.ArgumentParser", return_value=parser) as _ap, patch(
            "argparse.ArgumentParser.add_argument"
        ) as _aa:
            self.assertEqual(pim.setParser(), parser)
            _ap.assert_called_once_with(prog="prepareIni.py")
        parser.add_argument.assert_has_calls(
            [
                call(
                    "inifile",
                    metavar="inifilename",
                    help="the filename of the ini file where to write the configuration",
                ),
                call(
                    "outputfolder",
                    help="the name of the output folder that will contain the results",
                ),
                call(
                    "numodel",
                    choices=["3nu", "2nu", "3+1", "2+1", "1+1"],
                    help="define the neutrino model that must be used",
                ),
                call(
                    "--collint_diagonal_zero",
                    action="store_true",
                    help="set to zero the diagonal contributions of collision terms",
                ),
                call(
                    "--collint_offdiag_nodamp",
                    action="store_true",
                    help="use full integrals instead of damping terms for the off-diagonal collision terms",
                ),
                call(
                    "--collint_damping_type",
                    choices=["zero", "yyyw", "McKellar:1992ja"],
                    default="yyyw",
                    help="define the scheme for the off-diagonal contribution of collision integrals."
                    + "Use 'zero' to ignore all the off-diagonal components, "
                    + "'yyyw' to use the expressions from YYYW or "
                    + "'McKellar:1992ja' for the expressions from the paper McKellar:1992ja",
                ),
                call(
                    "--collint_d_no_nue",
                    action="store_true",
                    help="disable nue contributions to diagonal damping terms",
                ),
                call(
                    "--collint_d_no_nunu",
                    action="store_true",
                    help="disable nunu contributions to diagonal damping terms",
                ),
                call(
                    "--collint_od_no_nue",
                    action="store_true",
                    help="disable nue contributions to off-diagonal damping terms",
                ),
                call(
                    "--collint_od_no_nunu",
                    action="store_true",
                    help="disable nunu contributions to off-diagonal damping terms",
                ),
                call(
                    "--qed_corrections",
                    choices=["no", "o2", "o3", "o2ln", "o3ln"],
                    default="o3",
                    help="define which terms must be included "
                    + "for the finite-temperature QED corrections "
                    + "[O(e^2), O(e^2)+O(e^3) - default, O(e^2)+ln, O(e^2)+O(e^3)ln, or none]",
                ),
                call(
                    "--ordering",
                    choices=["NO", "IO"],
                    default="NO",
                    help="define the mass ordering for the three active neutrinos "
                    + "(not used if you explicitely give the mixing parameters, "
                    + "only if you select default values "
                    + "by Valencia, Bari or NuFit global fit)",
                ),
                call(
                    "--default_active",
                    choices=["Bari", "NuFit", "VLC", "None"],
                    default="VLC",
                    help="define the mixing parameters for the active neutrinos as obtained "
                    + "from the Valencia global fit "
                    + "(doi:10.1016/j.physletb.2018.06.019, default), "
                    + "Bari group (doi:10.1016/j.ppnp.2018.05.005) or "
                    + "NuFit analysis (doi:10.1007/JHEP01(2019)106)",
                ),
                call(
                    "--default_sterile",
                    choices=["Gariazzo&al", "None"],
                    default="Gariazzo&al",
                    help="define the active-sterile mixing parameters as obtained "
                    + "from the Gariazzo et al. global fit (with th24=th34=0)",
                ),
                call(
                    "--sinsq",
                    dest="use_sinsq",
                    action="store_false",
                    help=r"use the $\sin^2$ of the mixing angles as input",
                ),
                call(
                    "--dm21",
                    type=float,
                    default=0.0,
                    help=r"define $\Delta m^2_{21}$",
                ),
                call(
                    "--dm31",
                    type=float,
                    default=0.0,
                    help=r"define $\Delta m^2_{31}$ "
                    + "(pass negative value for inverted ordering)",
                ),
                call(
                    "--dm41",
                    type=float,
                    default=0.0,
                    help=r"define $\Delta m^2_{41}$",
                ),
                call(
                    "--th12",
                    type=float,
                    default=0.0,
                    help=r"define $\theta_{12}$ or $\sin^2 \theta_{12}$",
                ),
                call(
                    "--th13",
                    type=float,
                    default=0.0,
                    help=r"define $\theta_{13}$ or $\sin^2 \theta_{13}$",
                ),
                call(
                    "--th14",
                    type=float,
                    default=0.0,
                    help=r"define $\theta_{14}$ or $\sin^2 \theta_{14}$",
                ),
                call(
                    "--th23",
                    type=float,
                    default=0.0,
                    help=r"define $\theta_{23}$ or $\sin^2 \theta_{23}$",
                ),
                call(
                    "--th24",
                    type=float,
                    default=0.0,
                    help=r"define $\theta_{24}$ or $\sin^2 \theta_{24}$",
                ),
                call(
                    "--th34",
                    type=float,
                    default=0.0,
                    help=r"define $\theta_{34}$ or $\sin^2 \theta_{34}$",
                ),
                call(
                    "--Trh",
                    type=float,
                    default=25.0,
                    help=r"define $\T_{RH}$ (reheating temperature in MeV)",
                ),
                call(
                    "-V",
                    "--verbose",
                    type=int,
                    default=1,
                    help="define the verbosity of the code",
                ),
                call(
                    "--verbose_deriv_freq",
                    type=int,
                    default=100,
                    help="print a string stating the current position only after N derivatives",
                ),
                call(
                    "--Nx",
                    type=int,
                    default=200,
                    help="number of points to save in x",
                ),
                call("--x_in", type=float, default=0.01, help="initial value of x"),
                call("--x_fin", type=float, default=35, help="final value of x"),
                call("--Ny", type=int, default=30, help="number of total points in y"),
                call(
                    "--Nylog",
                    type=int,
                    default=5,
                    help="number of log-spaced points between y_in and y_cen",
                ),
                call("--y_min", type=float, default=0.01, help="minimum value of y"),
                call(
                    "--y_cen",
                    type=float,
                    default=1,
                    help="value of y where to switch between log- and linear spacing",
                ),
                call("--y_max", type=float, default=20, help="maximum value of y"),
                call(
                    "--dlsoda_atol",
                    type=float,
                    default=1e-6,
                    help="absolute tolerance for all the differential equations in DLSODA. "
                    + "See also dlsoda_atol_z, dlsoda_atol_d, dlsoda_atol_o",
                ),
                call(
                    "--dlsoda_atol_z",
                    type=float,
                    default=1e-6,
                    help="absolute tolerance for the dz/dx, dw/dx "
                    + "differential equations in DLSODA. "
                    + "See also dlsoda_atol, dlsoda_atol_d, dlsoda_atol_o",
                ),
                call(
                    "--dlsoda_atol_d",
                    type=float,
                    default=1e-6,
                    help="absolute tolerance for the differential equations drho_{ii}/dx"
                    + "of the diagonal matrix elements in DLSODA. "
                    + "See also dlsoda_atol, dlsoda_atol_z, dlsoda_atol_o",
                ),
                call(
                    "--dlsoda_atol_o",
                    type=float,
                    default=1e-6,
                    help="absolute tolerance for the differential equations drho_{ij}/dx"
                    + "of the off-diagonal matrix elements in DLSODA. "
                    + "See also dlsoda_atol, dlsoda_atol_z, dlsoda_atol_d",
                ),
                call(
                    "--dlsoda_rtol",
                    type=float,
                    default=1e-6,
                    help="relative tolerance for DLSODA",
                ),
                call(
                    "--save_BBN",
                    action="store_true",
                    help="enable saving the output for PArthENoPE",
                ),
                call(
                    "--save_energy_entropy",
                    action="store_true",
                    help="enable saving the evolution of the energy density "
                    + "and entropy for each component",
                ),
                call(
                    "--save_fd",
                    action="store_true",
                    help="enable saving the y grid and the corresponding Fermi-Dirac to fd.dat",
                ),
                call(
                    "--save_Neff",
                    action="store_true",
                    help="enable saving the evolution of Neff",
                ),
                call(
                    "--save_nuDens",
                    action="store_true",
                    help="enable saving the evolution of the full neutrino density matrix",
                ),
                call(
                    "--save_number",
                    action="store_true",
                    help="enable saving the evolution of the number density for each component",
                ),
                call(
                    "--save_z",
                    action="store_true",
                    help="enable saving the evolution of the photon temperature z",
                ),
                call(
                    "--no_GL",
                    action="store_true",
                    help="do not use the Gauss-Laguerre method for integrals "
                    + "and for spacing the y points",
                ),
            ],
            any_order=True,
        )

    def test_parsing(self):
        """test that there are no errors in the parsing process"""
        parser = pim.setParser()
        with self.assertRaises(SystemExit):
            args = parser.parse_args(["inifile", "outdir"])
        with self.assertRaises(SystemExit):
            args = parser.parse_args(["inifile", "outdir", "abc"])
        with self.assertRaises(SystemExit):
            args = parser.parse_args(["inifile", "outdir", "3nu", "abc"])
        for p in ["3nu", "2nu", "3+1", "2+1", "1+1"]:
            args = parser.parse_args(["inifile", "outdir", p])
        baseargs = ["inifile", "outdir", "3nu"]
        args = parser.parse_args(baseargs)
        self.assertEqual(args.ordering, "NO")
        self.assertEqual(args.default_active, "VLC")
        self.assertEqual(args.default_sterile, "Gariazzo&al")
        self.assertEqual(args.use_sinsq, True)
        self.assertEqual(args.dm21, 0.0)
        self.assertEqual(args.dm31, 0.0)
        self.assertEqual(args.dm41, 0.0)
        self.assertEqual(args.th12, 0.0)
        self.assertEqual(args.th13, 0.0)
        self.assertEqual(args.th14, 0.0)
        self.assertEqual(args.th23, 0.0)
        self.assertEqual(args.th24, 0.0)
        self.assertEqual(args.th34, 0.0)
        self.assertEqual(args.verbose, 1)
        self.assertEqual(args.verbose_deriv_freq, 100)
        self.assertEqual(args.Nx, 200)
        self.assertEqual(args.x_in, 0.01)
        self.assertEqual(args.x_fin, 35)
        self.assertEqual(args.Ny, 30)
        self.assertEqual(args.Nylog, 5)
        self.assertEqual(args.y_min, 0.01)
        self.assertEqual(args.y_cen, 1)
        self.assertEqual(args.y_max, 20)
        self.assertEqual(args.dlsoda_rtol, 1.0e-6)
        self.assertEqual(args.dlsoda_atol, 1.0e-6)
        self.assertEqual(args.dlsoda_atol_z, 1.0e-6)
        self.assertEqual(args.dlsoda_atol_d, 1.0e-6)
        self.assertEqual(args.dlsoda_atol_o, 1.0e-6)
        self.assertEqual(args.no_GL, False)
        self.assertEqual(args.collint_diagonal_zero, False)
        self.assertEqual(args.collint_offdiag_nodamp, False)
        self.assertEqual(args.collint_damping_type, "yyyw")
        self.assertEqual(args.collint_d_no_nue, False)
        self.assertEqual(args.collint_d_no_nunu, False)
        self.assertEqual(args.collint_od_no_nue, False)
        self.assertEqual(args.collint_od_no_nunu, False)
        self.assertEqual(args.save_BBN, False)
        self.assertEqual(args.save_energy_entropy, False)
        self.assertEqual(args.save_fd, False)
        self.assertEqual(args.save_Neff, False)
        self.assertEqual(args.save_nuDens, False)
        self.assertEqual(args.save_number, False)
        self.assertEqual(args.save_z, False)
        for l in [
            ["--ordering=NO"],
            ["--ordering=IO"],
            ["--default_active=None"],
            ["--default_active=Bari"],
            ["--default_active=NuFit"],
            ["--default_active=VLC"],
            ["--default_sterile=None"],
            ["--default_sterile=Gariazzo&al"],
            ["--sinsq"],
            ["--dm21=1.23"],
            ["--dm31=1.23"],
            ["--dm41=1.23"],
            ["--th12=0.01"],
            ["--th13=0.01"],
            ["--th14=0.01"],
            ["--th23=0.01"],
            ["--th24=0.01"],
            ["--th34=0.01"],
            ["-V=2"],
            ["--verbose=2"],
            ["--verbose_deriv_freq=200"],
            ["--Nx=100"],
            ["--x_in=0.1"],
            ["--x_fin=40.5"],
            ["--Ny=30"],
            ["--Nylog=5"],
            ["--y_min=0.02"],
            ["--y_cen=2.2"],
            ["--y_max=30.6"],
            ["--dlsoda_rtol=1.3e-5"],
            ["--dlsoda_atol=1.3e-4"],
            ["--dlsoda_atol_z=1.2e-6"],
            ["--dlsoda_atol_d=1.2e-6"],
            ["--dlsoda_atol_o=1.2e-6"],
            ["--no_GL"],
            ["--save_BBN"],
            ["--save_energy_entropy"],
            ["--save_fd"],
            ["--save_Neff"],
            ["--save_nuDens"],
            ["--save_number"],
            ["--save_z"],
        ]:
            args = parser.parse_args(baseargs + l)
        for l in [
            ["--someoption"],
            ["--ordering=abc"],
            ["--default_active=abc"],
            ["--dm21=abc"],
            ["--dm31=abc"],
            ["--dm41=abc"],
            ["--th12=abc"],
            ["--th13=abc"],
            ["--th14=abc"],
            ["--th23=abc"],
            ["--th24=abc"],
            ["--th34=abc"],
            ["--verbose=abc"],
            ["--verbose_deriv_freq=abc"],
            ["--Nx=abc"],
            ["--x_in=abc"],
            ["--x_fin=abc"],
            ["--Ny=abc"],
            ["--Nylog=abc"],
            ["--y_min=abc"],
            ["--y_cen=abc"],
            ["--y_max=abc"],
            ["--dlsoda_rtol=abc"],
            ["--dlsoda_atol=abc"],
            ["--dlsoda_atol_z=abc"],
            ["--dlsoda_atol_d=abc"],
            ["--dlsoda_atol_o=abc"],
            ["--no_GL=a"],
            ["--save_BBN=a"],
            ["--save_energy_entropy=a"],
            ["--save_fd=a"],
            ["--save_Neff=a"],
            ["--save_nuDens=a"],
            ["--save_number=a"],
            ["--save_z=a"],
        ]:
            with self.assertRaises(SystemExit):
                args = parser.parse_args(baseargs + l)
        values = pim.getIniValues(args)
        pim.writeIni(self.testIni, values)

    def test_oscParams(self):
        """test oscParams"""
        args = Namespace(
            **{
                "use_sinsq": True,
                "dm21": 1.23,
                "th12": 0.01,
            }
        )
        for a in ["a+s", "as", "1p1", "1+1"]:
            args.numodel = a
            self.assertEqual(
                pim.oscParams(args),
                {
                    "use_sinsq": "T",
                    "nnu": 2,
                    "factors": [1, 1],
                    "sterile": [False, True],
                    "dm41": 0.0,
                    "th14": 0.0,
                    "th24": 0.0,
                    "th34": 0.0,
                    "dm31": 0.0,
                    "th13": 0.0,
                    "th23": 0.0,
                    "dm21": 1.23,
                    "th12": 0.01,
                },
            )
        for a in ["2+0", "2nu", "2"]:
            args.numodel = a
            self.assertEqual(
                pim.oscParams(args),
                {
                    "use_sinsq": "T",
                    "nnu": 2,
                    "factors": [1, 2],
                    "sterile": [False, False],
                    "dm41": 0.0,
                    "th14": 0.0,
                    "th24": 0.0,
                    "th34": 0.0,
                    "dm31": 0.0,
                    "th13": 0.0,
                    "th23": 0.0,
                    "dm21": 1.23,
                    "th12": 0.01,
                },
            )
        args = Namespace(
            **{
                "use_sinsq": True,
                "dm21": 1.23,
                "th12": 0.01,
                "dm31": 0.0025,
                "th13": 0.02,
                "th23": 0.5,
            }
        )
        for a in ["2p1", "2+1"]:
            args.numodel = a
            self.assertEqual(
                pim.oscParams(args),
                {
                    "use_sinsq": "T",
                    "nnu": 3,
                    "factors": [1, 2, 1],
                    "sterile": [False, False, True],
                    "dm41": 0.0,
                    "th14": 0.0,
                    "th24": 0.0,
                    "th34": 0.0,
                    "dm31": 0.0025,
                    "th13": 0.02,
                    "th23": 0.5,
                    "dm21": 1.23,
                    "th12": 0.01,
                },
            )
        args = Namespace(
            **{
                "use_sinsq": False,
                "dm21": 1.23,
                "th12": 0.01,
                "dm31": 0.0025,
                "th13": 0.02,
                "th23": 0.5,
            }
        )
        for a in ["3p0", "3+0", "3nu", "3"]:
            args.numodel = a
            args.default_active = "g"
            args.ordering = "a"
            self.assertEqual(
                pim.oscParams(args),
                {
                    "use_sinsq": "F",
                    "nnu": 3,
                    "factors": [1, 1, 1],
                    "sterile": [False, False, False],
                    "dm41": 0.0,
                    "th14": 0.0,
                    "th24": 0.0,
                    "th34": 0.0,
                    "dm31": 0.0025,
                    "th13": 0.02,
                    "th23": 0.5,
                    "dm21": 1.23,
                    "th12": 0.01,
                },
            )
            for g in ["VLC", "Bari", "NuFit"]:
                for o in ["NO", "IO"]:
                    args.default_active = g
                    args.ordering = o
                    res = pim.default_osc_act[g][o.lower()].copy()
                    res.update(
                        {
                            "use_sinsq": "F",
                            "nnu": 3,
                            "factors": [1, 1, 1],
                            "sterile": [False, False, False],
                            "dm41": 0.0,
                            "th14": 0.0,
                            "th24": 0.0,
                            "th34": 0.0,
                        }
                    )
                    self.assertEqual(pim.oscParams(args), res)
        args = Namespace(
            **{
                "use_sinsq": True,
                "dm21": 8e-5,
                "th12": 0.01,
                "dm31": 0.0025,
                "th13": 0.02,
                "th23": 0.5,
                "dm41": 1.23,
                "th14": 0.1,
                "th24": 0.2,
                "th34": 0.3,
            }
        )
        for a in ["3p1", "3+1"]:
            args.numodel = a
            args.default_sterile = "g"
            args.default_active = "g"
            args.ordering = "a"
            self.assertEqual(
                pim.oscParams(args),
                {
                    "use_sinsq": "T",
                    "nnu": 4,
                    "factors": [1, 1, 1, 1],
                    "sterile": [False, False, False, True],
                    "dm41": 1.23,
                    "th14": 0.1,
                    "th24": 0.2,
                    "th34": 0.3,
                    "dm31": 0.0025,
                    "th13": 0.02,
                    "th23": 0.5,
                    "dm21": 8e-5,
                    "th12": 0.01,
                },
            )
            for g in ["Gariazzo&al"]:
                args.default_sterile = g
                args.default_active = "g"
                args.ordering = "a"
                res = pim.default_osc_ster[g].copy()
                res.update(
                    {
                        "use_sinsq": "T",
                        "nnu": 4,
                        "factors": [1, 1, 1, 1],
                        "sterile": [False, False, False, True],
                        "dm31": 0.0025,
                        "th13": 0.02,
                        "th23": 0.5,
                        "dm21": 8e-5,
                        "th12": 0.01,
                    }
                )
                self.assertEqual(pim.oscParams(args), res)
            args.default_sterile = "g"
            for g in ["VLC", "Bari", "NuFit"]:
                for o in ["NO", "IO"]:
                    args.default_active = g
                    args.ordering = o
                    res = pim.default_osc_act[g][o.lower()].copy()
                    res.update(
                        {
                            "use_sinsq": "T",
                            "nnu": 4,
                            "factors": [1, 1, 1, 1],
                            "sterile": [False, False, False, True],
                            "dm41": 1.23,
                            "th14": 0.1,
                            "th24": 0.2,
                            "th34": 0.3,
                        }
                    )
                    self.assertEqual(pim.oscParams(args), res)

    def test_getIniValues(self):
        """test getIniValues"""
        self.maxDiff = None
        args = Namespace(
            **{
                "verbose": "vb",
                "collint_diagonal_zero": False,
                "collint_offdiag_nodamp": False,
                "collint_damping_type": "yyyw",
                "Nx": 200,
                "x_in": 0.001,
                "x_fin": 35,
                "Trh": 123.0,
                "Ny": 24,
                "Nylog": 4,
                "y_min": 0.01,
                "y_cen": 1,
                "y_max": 20,
                "dlsoda_atol_z": 1e-6,
                "dlsoda_atol_d": 1e-6,
                "dlsoda_atol_o": 1e-7,
                "dlsoda_rtol": 1e-4,
                "outputfolder": "abcd",
                "verbose_deriv_freq": 123,
                "no_GL": True,
                "collint_d_no_nue": True,
                "collint_d_no_nunu": True,
                "collint_od_no_nue": True,
                "collint_od_no_nunu": True,
                "save_BBN": True,
                "save_energy_entropy": True,
                "save_fd": True,
                "save_Neff": True,
                "save_nuDens": True,
                "save_number": True,
                "save_z": True,
                "qed_corrections": "no",
            }
        )
        with patch(
            "prepareIni.oscParams",
            return_value={"factors": [2, 1], "sterile": [False, True]},
        ) as _op:
            values = pim.getIniValues(args)
        self.assertEqual(
            values,
            {
                "verbose": "vb",
                "factors": "nuFactor1 = %f\nnuFactor2 = %f" % (2, 1),
                "sterile": "sterile1 = F\nsterile2 = T",
                "collint_diagonal_zero": "F",
                "collint_offdiag_damping": "T",
                "collint_damping_type": 1,
                "collint_d_no_nue": "T",
                "collint_d_no_nunu": "T",
                "collint_od_no_nue": "T",
                "collint_od_no_nunu": "T",
                "Nx": 200,
                "x_in": 0.001,
                "x_fin": 35,
                "Trh": 123.0,
                "Ny": 24,
                "Nylog": 4,
                "y_min": 0.01,
                "y_cen": 1,
                "y_max": 20,
                "dlsoda_atol": "dlsoda_atol_z = %s\n" % 1e-6
                + "dlsoda_atol_d = %s\n" % 1e-6
                + "dlsoda_atol_o = %s\n" % 1e-7,
                "dlsoda_rtol": 1e-4,
                "folder": "abcd",
                "Nprintderivs": 123,
                "use_GL": "F",
                "save_BBN": "T",
                "save_energy_entropy": "T",
                "save_fd": "T",
                "save_Neff": "T",
                "save_nuDens": "T",
                "save_number": "T",
                "save_z": "T",
                "ftqed_temperature_corr": "F",
                "ftqed_ord3": "F",
                "ftqed_log_term": "F",
            },
        )
        args = Namespace(
            **{
                "verbose": "vb",
                "collint_diagonal_zero": True,
                "collint_offdiag_nodamp": True,
                "collint_damping_type": "zero",
                "Nx": 200,
                "x_in": 0.001,
                "x_fin": 35,
                "Trh": 123.0,
                "Ny": 24,
                "Nylog": 4,
                "y_min": 0.01,
                "y_cen": 1,
                "y_max": 20,
                "dlsoda_atol": 1e-5,
                "dlsoda_atol_z": 1e-6,
                "dlsoda_atol_d": 1e-6,
                "dlsoda_atol_o": 1e-6,
                "dlsoda_rtol": 1e-4,
                "collint_d_no_nue": False,
                "collint_d_no_nunu": False,
                "collint_od_no_nue": False,
                "collint_od_no_nunu": False,
                "outputfolder": "abcd",
                "verbose_deriv_freq": 123,
                "no_GL": False,
                "save_BBN": False,
                "save_energy_entropy": False,
                "save_fd": False,
                "save_Neff": False,
                "save_nuDens": False,
                "save_number": False,
                "save_z": False,
                "qed_corrections": "o2",
            }
        )
        with patch(
            "prepareIni.oscParams",
            return_value={"factors": [2, 1], "sterile": [False, True]},
        ) as _op:
            values = pim.getIniValues(args)
        self.assertEqual(
            values,
            {
                "verbose": "vb",
                "factors": "nuFactor1 = %f\nnuFactor2 = %f" % (2, 1),
                "sterile": "sterile1 = F\nsterile2 = T",
                "collint_diagonal_zero": "T",
                "collint_offdiag_damping": "F",
                "collint_damping_type": 0,
                "Nx": 200,
                "x_in": 0.001,
                "x_fin": 35,
                "Trh": 123.0,
                "Ny": 24,
                "Nylog": 4,
                "y_min": 0.01,
                "y_cen": 1,
                "y_max": 20,
                "collint_d_no_nue": "F",
                "collint_d_no_nunu": "F",
                "collint_od_no_nue": "F",
                "collint_od_no_nunu": "F",
                "dlsoda_atol": "dlsoda_atol_z = %s\n" % 1e-5
                + "dlsoda_atol_d = %s\n" % 1e-5
                + "dlsoda_atol_o = %s\n" % 1e-5,
                "dlsoda_rtol": 1e-4,
                "folder": "abcd",
                "Nprintderivs": 123,
                "use_GL": "T",
                "save_BBN": "F",
                "save_energy_entropy": "F",
                "save_fd": "F",
                "save_Neff": "F",
                "save_nuDens": "F",
                "save_number": "F",
                "save_z": "F",
                "ftqed_temperature_corr": "T",
                "ftqed_ord3": "F",
                "ftqed_log_term": "F",
            },
        )
        args.collint_damping_type = "McKellar:1992ja"
        args.qed_corrections = "o3"
        args.dlsoda_atol_z = 1e-7
        args.dlsoda_atol_d = 1e-6
        args.dlsoda_atol_o = 1e-6
        with patch(
            "prepareIni.oscParams",
            return_value={"factors": [2, 1], "sterile": [False, True]},
        ) as _op:
            values = pim.getIniValues(args)
        self.assertEqual(values["collint_damping_type"], 2)
        self.assertEqual(values["ftqed_temperature_corr"], "T")
        self.assertEqual(values["ftqed_ord3"], "T")
        self.assertEqual(values["ftqed_log_term"], "F")
        self.assertEqual(
            values["dlsoda_atol"],
            "dlsoda_atol_z = %s\n" % 1e-7
            + "dlsoda_atol_d = %s\n" % 1e-6
            + "dlsoda_atol_o = %s\n" % 1e-6,
        )
        args.collint_damping_type = "abcd"
        args.qed_corrections = "o2ln"
        args.dlsoda_atol_z = 1e-6
        args.dlsoda_atol_d = 1e-7
        args.dlsoda_atol_o = 1e-6
        with patch(
            "prepareIni.oscParams",
            return_value={"factors": [2, 1], "sterile": [False, True]},
        ) as _op:
            values = pim.getIniValues(args)
        self.assertEqual(values["collint_damping_type"], 1)
        self.assertEqual(values["ftqed_temperature_corr"], "T")
        self.assertEqual(values["ftqed_ord3"], "F")
        self.assertEqual(values["ftqed_log_term"], "T")
        self.assertEqual(
            values["dlsoda_atol"],
            "dlsoda_atol_z = %s\n" % 1e-6
            + "dlsoda_atol_d = %s\n" % 1e-7
            + "dlsoda_atol_o = %s\n" % 1e-6,
        )
        args.qed_corrections = "o3ln"
        with patch(
            "prepareIni.oscParams",
            return_value={"factors": [2, 1], "sterile": [False, True]},
        ) as _op:
            values = pim.getIniValues(args)
        self.assertEqual(values["ftqed_temperature_corr"], "T")
        self.assertEqual(values["ftqed_ord3"], "T")
        self.assertEqual(values["ftqed_log_term"], "T")

    def test_writeIni(self):
        """test writeIni"""
        values = {
            "verbose": "vb",
            "use_sinsq": "T",
            "nnu": 4,
            "dm41": 1.23,
            "th14": 0.1,
            "th24": 0.2,
            "th34": 0.3,
            "dm31": 0.0025,
            "th13": 0.02,
            "th23": 0.5,
            "dm21": 8e-5,
            "th12": 0.01,
            "Trh": 123.0,
            "factors": "nuFactor1 = %f\nnuFactor2 = %f" % (2, 1),
            "sterile": "sterile1 = F\nsterile2 = T",
            "collint_diagonal_zero": "F",
            "collint_offdiag_damping": "T",
            "collint_damping_type": "1",
            "collint_d_no_nue": "F",
            "collint_d_no_nunu": "F",
            "collint_od_no_nue": "F",
            "collint_od_no_nunu": "F",
            "Nx": 200,
            "x_in": 0.001,
            "x_fin": 35,
            "Ny": 24,
            "Nylog": 4,
            "y_min": 0.01,
            "y_cen": 1,
            "y_max": 20,
            "ftqed_temperature_corr": "F",
            "ftqed_ord3": "T",
            "ftqed_log_term": "T",
            "dlsoda_atol": "dlsoda_atol_z = %s\n" % 1e-5
            + "dlsoda_atol_d = %s\n" % 1e-5
            + "dlsoda_atol_o = %s\n" % 1e-5,
            "dlsoda_rtol": 1e-4,
            "folder": "abcd/",
            "Nprintderivs": 123,
            "use_GL": "T",
            "save_BBN": "F",
            "save_energy_entropy": "F",
            "save_fd": "F",
            "save_Neff": "F",
            "save_nuDens": "F",
            "save_number": "F",
            "save_z": "F",
        }
        with self.assertRaises(OSError):
            pim.writeIni("/nonexistent/cant/write/this.ini", values)
        pim.writeIni(self.testIni, values)
        with open(self.testIni) as _f:
            lines = _f.readlines()
        for l in [
            "flavorNumber = 4",
            "givesinsq = T",
            "theta12= 0.01",
            "dm21 = %s" % 8e-5,
            "theta13 = 0.02",
            "theta23 = 0.5",
            "dm31 = 0.0025",
            "theta14 = 0.1",
            "theta24 = 0.2",
            "theta34 = 0.3",
            "dm41 = 1.23",
            "collint_diagonal_zero = F",
            "collint_offdiag_damping = T",
            "collint_damping_type = 1",
            "ftqed_temperature_corr = F",
            "collint_d_no_nue = F",
            "collint_d_no_nunu = F",
            "collint_od_no_nue = F",
            "collint_od_no_nunu = F",
            "ftqed_ord3 = T",
            "ftqed_log_term = T",
            "Nx = 200",
            "x_in = 0.001",
            "x_fin = 35",
            "use_gauss_laguerre = T",
            "Ny = 24",
            "Nylog = 4",
            "y_min = 0.01",
            "y_cen = 1",
            "y_max = 20",
            "outputFolder = abcd/",
            "checkpoint = T",
            "save_BBN = F",
            "save_fd = F",
            "save_Neff = F",
            "save_nuDens_evolution = F",
            "save_z_evolution = F",
            "save_energy_entropy_evolution = F",
            "save_number_evolution = F",
            "dlsoda_rtol = %s" % 1e-4,
            "verbose = vb",
            "Nprintderivs = 123",
        ]:
            self.assertIn(l + "\n", lines)
        with open(self.testIni) as _f:
            content = _f.read()
        for c in [
            "nuFactor1 = %f\nnuFactor2 = %f" % (2, 1),
            "sterile1 = F\nsterile2 = T",
            "dlsoda_atol_z = %s\n" % 1e-5
            + "dlsoda_atol_d = %s\n" % 1e-5
            + "dlsoda_atol_o = %s\n" % 1e-5,
        ]:
            self.assertIn(c, content)


if __name__ == "__main__":
    unittest.main()
