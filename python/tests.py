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


class TestFortepianoOutput(unittest.TestCase):
    """Testing some content of the FortepianoOutput module"""

    def test_attributes(self):
        """test attributes"""
        self.assertIsInstance(fpom.colors, list)
        self.assertIsInstance(fpom.styles, list)
        self.assertIsInstance(fpom.markers, list)
        self.assertIsInstance(fpom.PISQD15, float)

    def test_finalizePlot(self):
        """test finalizePlot"""

    def test_stripRepeated(self):
        """test stripRepeated"""


class TestFortEPiaNORun(unittest.TestCase):
    """Testing the FortEPiaNORun class"""

    @classmethod
    def setUpClass(self):
        """Set maxDiff to None and load the output of explanatory.ini"""
        self.maxDiff = None
        self.explanatory = fpom.FortEPiaNORun("output/", label="label")

    def assertEqualArray(self, a, b):
        # type: (np.ndarray, np.ndarray) -> None
        """Assert that two np.ndarrays (a, b) are equal
        using np.all(a == b)

        Parameters:
            a, b: the two np.ndarrays to compare
        """
        self.assertTrue(np.allclose(a, b, equal_nan=True))

    def runAllPlots(self, run):
        """Call all the plotting function from FortEPiaNORun"""
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
        self.assertFalse(hasattr(run, "eldens"))
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
        self.assertTrue(hasattr(fpom.FortEPiaNORun, "interpolateRhoIJ"))

    def test_interpolateRhoIJ_x(self):
        """test interpolateRhoIJ_x"""
        self.assertTrue(hasattr(fpom.FortEPiaNORun, "interpolateRhoIJ_x"))

    def test_printTableLine(self):
        """test printTableLine"""
        self.assertTrue(hasattr(fpom.FortEPiaNORun, "printTableLine"))

    def test_plotFD(self):
        """test plotFD"""
        self.assertTrue(hasattr(fpom.FortEPiaNORun, "plotFD"))

    def test_plotZ(self):
        """test plotZ"""
        self.assertTrue(hasattr(fpom.FortEPiaNORun, "plotZ"))

    def test_plotW(self):
        """test plotW"""
        self.assertTrue(hasattr(fpom.FortEPiaNORun, "plotW"))

    def test_plotZoverW(self):
        """test plotZoverW"""
        self.assertTrue(hasattr(fpom.FortEPiaNORun, "plotZoverW"))

    def test_plotDeltaZ(self):
        """test plotDeltaZ"""
        self.assertTrue(hasattr(fpom.FortEPiaNORun, "plotDeltaZ"))

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
        self.assertTrue(hasattr(fpom.FortEPiaNORun, "doAllPlots"))

    def test_integrateRho_yn(self):
        """test integrateRho_yn"""
        self.assertTrue(hasattr(fpom.FortEPiaNORun, "integrateRho_yn"))


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
