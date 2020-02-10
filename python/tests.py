#!/usr/bin/python
import os
import shutil
import sys
import numpy as np
import matplotlib

matplotlib.use("agg")
import matplotlib.pyplot as plt

if sys.version_info[0] < 3:
    import unittest2 as unittest
    from mock import call, patch, MagicMock

    USE_AUTOSPEC_CLASS = False
else:
    import unittest
    from unittest.mock import call, patch, MagicMock

    USE_AUTOSPEC_CLASS = True


import dogrid3p1 as dgm
import fortepianoOutput as fpom
import prepareIni as pim

try:
    FileNotFoundError
except NameError:
    FileNotFoundError = IOError

allPlots = True


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


class TestDogrid3p1(unittest.TestCase):
    """Testing the dogrid3p1 module"""

    def test_sortFiles(self):
        """test sortFiles"""

    def test_Namespace(self):
        """test Namespace"""

    def test_safegetattr(self):
        """test savegetattr"""

    def test_setParser(self):
        """test setParser"""

    def test_write_grid_cfg(self):
        """test write_grid_cfg"""

    def test_read_grid_cfg(self):
        """test read_grid_cfg"""

    def test_contourplot(self):
        """test contourplot"""
        self.assertTrue(hasattr(dgm, "contourplot"))

    def test_call_fill(self):
        """test call_fill"""
        self.assertTrue(hasattr(dgm, "call_fill"))

    def test_call_plot(self):
        """test call_plot"""
        self.assertTrue(hasattr(dgm, "call_plot"))

    def test_call_prepare(self):
        """test call_prepare"""
        self.assertTrue(hasattr(dgm, "call_prepare"))

    def test_call_read(self):
        """test call_read"""
        self.assertTrue(hasattr(dgm, "call_read"))

    def test_call_run(self):
        """test call_run"""
        self.assertTrue(hasattr(dgm, "call_run"))

    def test_call_ternary(self):
        """test call_ternary"""
        self.assertTrue(hasattr(dgm, "call_ternary"))


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
            _tig.assert_called_once_with(rect=(-0.035, -0.04, 1.025, 1.04))
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
                "fname", legcol=3, lloc="upper right",
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

    def test_example(self):
        """test an example with FortEPiaNORun from explanatory.ini"""
        folder = "output/"
        run = fpom.FortEPiaNORun(folder, label="label")
        fc = np.loadtxt("%s/fd.dat" % folder)
        self.assertEqualArray(run.yv, fc[:, 0])
        self.assertEqualArray(run.fd, fc[:, 1])
        fc = np.loadtxt("%s/z.dat" % folder)
        self.assertEqualArray(run.zdat, fc)
        fc = np.loadtxt("%s/Neff.dat" % folder)
        self.assertEqualArray(run.Neffdat, fc)
        fc = np.loadtxt("%s/energyDensity.dat" % folder)
        self.assertEqualArray(run.endens, fc)
        fc = np.loadtxt("%s/entropy.dat" % folder)
        self.assertEqualArray(run.entropy, fc)
        self.assertTrue(hasattr(run, "resume"))
        self.assertTrue(run.hasResume)
        self.assertTrue(np.isclose(run.Neff, 3.0430, atol=1e-4))
        self.assertTrue(np.isclose(run.wfin, 1.09659, atol=1e-5))
        self.assertTrue(np.isclose(run.zfin, 1.53574, atol=1e-5))
        self.assertIsInstance(run.deltarhofin, list)
        self.assertTrue(np.isclose(run.deltarhofin[0], 0.4667, atol=1e-4))
        self.assertTrue(np.isclose(run.deltarhofin[1], 0.4639, atol=1e-4))
        self.assertTrue(np.isclose(run.deltarhofin[2], 0.4628, atol=1e-4))
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
        # now just do plots in order to see that everything works till the end
        self.runAllPlots(run)

    def test_failing(self):
        """test few failing examples with FortEPiaNORun"""
        folder = "output/no/"
        if os.path.exists(folder):
            shutil.rmtree(folder)
        run = fpom.FortEPiaNORun(folder)
        self.assertFalse(hasattr(run, "yv"))
        self.assertFalse(hasattr(run, "fd"))
        self.assertFalse(hasattr(run, "zdat"))
        self.assertFalse(hasattr(run, "Neffdat"))
        self.assertFalse(hasattr(run, "endens"))
        self.assertFalse(hasattr(run, "entropy"))
        self.assertFalse(hasattr(run, "Neff"))
        self.assertFalse(hasattr(run, "wfin"))
        self.assertFalse(hasattr(run, "zfin"))
        self.assertFalse(hasattr(run, "rho"))
        self.assertFalse(hasattr(run, "rhoM"))
        self.assertFalse(hasattr(run, "resume"))
        self.assertFalse(hasattr(run, "hasResume"))
        self.assertFalse(hasattr(run, "deltarhofin"))
        self.runAllPlots(run)

        # repeat creating some bad resume file, e.g. with nans
        os.makedirs(folder)
        with open("%s/resume.dat" % folder, "w") as _f:
            _f.write("final w =  NaN\nfinal z =  1.5\nNeff    =  \n")
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
            resume = _f.readlines()
        self.assertEqual(run.resume, resume)
        self.assertTrue(run.hasResume)
        self.assertEqualArray(run.deltarhofin, [np.nan, np.nan, np.nan])
        self.runAllPlots(run)

    def test_init(self):
        """test __init__"""
        folder = "output/"
        with patch("fortepianoOutput.FortEPiaNORun.doAllPlots") as _pl, patch(
            "fortepianoOutput.FortEPiaNORun.printTableLine"
        ) as _ptl:
            run = fpom.FortEPiaNORun(folder)
            self.assertEqual(_pl.call_count, 0)
            _ptl.assert_called_once()
        self.assertEqual(run.folder, folder)
        self.assertTrue(run.full)
        self.assertEqual(run.label, "")
        self.assertEqual(run.nnu, 3)
        self.assertTrue(run.verbose)

        with patch("fortepianoOutput.FortEPiaNORun.doAllPlots") as _pl, patch(
            "fortepianoOutput.FortEPiaNORun.printTableLine"
        ) as _ptl:
            run = fpom.FortEPiaNORun(
                folder, label="l", full=False, nnu=2, verbose=False, plots=True
            )
            _pl.assert_called_once()
            _ptl.assert_called_once()
        self.assertEqual(run.folder, folder)
        self.assertFalse(run.full)
        self.assertEqual(run.label, "l")
        self.assertEqual(run.nnu, 2)
        self.assertFalse(run.verbose)

    def test_interpolateRhoIJ(self):
        """test interpolateRhoIJ"""
        run = fpom.FortEPiaNORun("output/nonexistent")
        run.rho = np.asarray([])
        run.rhoM = np.asarray([])
        with self.assertRaises(IndexError):
            run.interpolateRhoIJ(0, 0, 5)
        run.rho = np.asarray(
            [[[None, None], [None, None],], [[None, None], [None, None],],],
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
            [[[None, None], [None, None],], [[None, None], [None, None],],],
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
            run.interpolateRhoIJ_x(0, 0, 2, y2=True), [[0, 1, 2], [0, 22, 4 * 32]],
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
        self.assertTrue(hasattr(fpom.FortEPiaNORun, "printTableLine"))

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
                _p.call_args[1], {"label": "label", "ls": "-", "marker": ".", "c": "k",}
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
                _p.call_args[1], {"label": "l", "ls": ":", "marker": ".", "c": "r",}
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
            self.assertEqual(_p.call_args[1], {"label": "label", "ls": "-", "c": "k",})
            _xl.assert_called_once_with("$x$")
            _yl.assert_called_once_with(r"$z$")
            _xs.assert_called_once_with("log")
            self.explanatory.plotZ(ls=":", lc="r", lab="l")
            self.assertEqual(_p.call_args[1], {"label": "l", "ls": ":", "c": "r",})
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
            self.assertEqual(_p.call_args[1], {"label": "label", "ls": "-", "c": "k",})
            _xl.assert_called_once_with("$x$")
            _yl.assert_called_once_with(r"$w$")
            _xs.assert_called_once_with("log")
            self.explanatory.plotW(ls=":", lc="r", lab="l")
            self.assertEqual(_p.call_args[1], {"label": "l", "ls": ":", "c": "r",})
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
            self.assertEqual(_p.call_args[1], {"label": "label", "ls": "-", "c": "k",})
            _xl.assert_called_once_with("$x$")
            _yl.assert_called_once_with(r"$z/w$")
            _xs.assert_called_once_with("log")
            self.explanatory.plotZoverW(ls=":", lc="r", lab="l")
            self.assertEqual(_p.call_args[1], {"label": "l", "ls": ":", "c": "r",})
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
            self.assertEqual(_p.call_args[1], {"label": "label", "ls": "-", "c": "k",})
            _xl.assert_called_once_with("$x$")
            _yl.assert_called_once_with(r"$z-z_{\rm ref}$")
            _xs.assert_called_once_with("log")
            run.plotDeltaZ(run1)
            self.assertEqualArray(_p.call_args[0], np.asarray([xv, -0.1 * yv]))
            run1.plotDeltaZ(run)
            self.assertEqualArray(_p.call_args[0], np.asarray([xv, 0.1 * yv]))
            self.explanatory.plotDeltaZ(run, lab="l", ls=":", lc="r")
            self.assertEqual(_p.call_args[1], {"label": "l", "ls": ":", "c": "r",})
            self.failedRun.plotDeltaZ(run)
            self.assertEqualArray(_p.call_args[0], [[np.nan, np.nan, np.nan]] * 2)
            run.plotDeltaZ(self.failedRun)
            self.assertEqualArray(_p.call_args[0], [xv, [np.nan] * len(xv)])

    def test_plotRhoDiag(self):
        """test plotRhoDiag"""
        self.assertTrue(hasattr(fpom.FortEPiaNORun, "plotRhoDiag"))

    def test_plotdRhoDiag(self):
        """test plotdRhoDiag"""
        self.assertTrue(hasattr(fpom.FortEPiaNORun, "plotdRhoDiag"))

    def test_plotRhoOffDiag(self):
        """test plotRhoOffDiag"""
        self.assertTrue(hasattr(fpom.FortEPiaNORun, "plotRhoOffDiag"))

    def test_plotdRhoOffDiag(self):
        """test plotdRhoOffDiag"""
        self.assertTrue(hasattr(fpom.FortEPiaNORun, "plotdRhoOffDiag"))

    def test_plotRhoFin(self):
        """test plotRhoFin"""
        self.assertTrue(hasattr(fpom.FortEPiaNORun, "plotRhoFin"))

    def test_plotRhoX(self):
        """test plotRhoX"""
        self.assertTrue(hasattr(fpom.FortEPiaNORun, "plotRhoX"))

    def test_plotRhoDiagY(self):
        """test plotRhoDiagY"""
        self.assertTrue(hasattr(fpom.FortEPiaNORun, "plotRhoDiagY"))

    def test_plotdRhoDiagY(self):
        """test plotdRhoDiagY"""
        self.assertTrue(hasattr(fpom.FortEPiaNORun, "plotdRhoDiagY"))

    def test_plotRhoOffDiagY(self):
        """test plotRhoOffDiagY"""
        self.assertTrue(hasattr(fpom.FortEPiaNORun, "plotRhoOffDiagY"))

    def test_plotdRhoOffDiagY(self):
        """test plotdRhoOffDiagY"""
        self.assertTrue(hasattr(fpom.FortEPiaNORun, "plotdRhoOffDiagY"))

    def test_plotNeff(self):
        """test plotNeff"""
        self.assertTrue(hasattr(fpom.FortEPiaNORun, "plotNeff"))

    def test_plotEnergyDensity(self):
        """test plotEnergyDensity"""
        self.assertTrue(hasattr(fpom.FortEPiaNORun, "plotEnergyDensity"))

    def test_plotEntropy(self):
        """test plotEntropy"""
        self.assertTrue(hasattr(fpom.FortEPiaNORun, "plotEntropy"))

    def test_doAllPlots(self):
        """test doAllPlots"""
        run = self.explanatory
        with patch("matplotlib.pyplot.close") as _cl, patch(
            "fortepianoOutput.finalizePlot"
        ) as _fp, patch("fortepianoOutput.FortEPiaNORun.plotZ") as _plZ, patch(
            "fortepianoOutput.FortEPiaNORun.plotW"
        ) as _plW, patch(
            "fortepianoOutput.FortEPiaNORun.plotRhoDiagY"
        ) as _plRdy, patch(
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
            self.assertEqual(_plRdy.call_count, 6)
            for i in range(run.nnu):
                _plRdy.assert_any_call(i, 5.0, fpom.styles[i], lc=fpom.colors[i])
                _plRdy.assert_any_call(
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
                        xlab="$x$",
                        ylab=r"$\rho$",
                        xscale="log",
                        yscale="log",
                    ),
                    call(
                        "%s/rho_mass_diag.pdf" % run.folder,
                        xlab="$x$",
                        ylab=r"$\rho$",
                        xscale="log",
                        yscale="log",
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
                ]
            )

        with patch("matplotlib.pyplot.close") as _cl, patch(
            "fortepianoOutput.finalizePlot"
        ) as _fp, patch("fortepianoOutput.FortEPiaNORun.plotZ") as _plZ, patch(
            "fortepianoOutput.FortEPiaNORun.plotW"
        ) as _plW, patch(
            "fortepianoOutput.FortEPiaNORun.plotRhoDiagY"
        ) as _plRdy, patch(
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
            self.assertEqual(_plRdy.call_count, 6)
            for i in range(run.nnu):
                _plRdy.assert_any_call(i, 2.5, fpom.styles[i], lc=fpom.colors[i])
                _plRdy.assert_any_call(
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
                        xlab="$x$",
                        ylab=r"$\rho$",
                        xscale="log",
                        yscale="log",
                    ),
                    call(
                        "%s/rho_mass_diag.pdf" % run.folder,
                        xlab="$x$",
                        ylab=r"$\rho$",
                        xscale="log",
                        yscale="log",
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
            [[[None, None], [None, None],], [[None, None], [None, None],],],
        )
        run.rhoM = np.asarray(
            [[[None, None], [None, None],], [[None, None], [None, None],],],
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

    def test_setParser(self):
        """test setParser"""

    def test_oscParams(self):
        """test oscParams"""

    def test_getIniValues(self):
        """test getIniValues"""

    def test_writeIni(self):
        """test writeIni"""


if __name__ == "__main__":
    unittest.main()
