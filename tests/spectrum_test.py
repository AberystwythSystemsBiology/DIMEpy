import unittest
from dimepy import Spectrum
import numpy as np
import os

scrpt_dir = os.path.dirname(os.path.abspath(__file__))

mzml_fp = os.path.join(scrpt_dir, "data/example.mzML")

class SpectrumTest(unittest.TestCase):

    def test_loading(self):
        spectrum = Spectrum(mzml_fp, "test")
        spectrum.load_scans()

    def test_negative_polarity(self):
        # The entire spectrum is positive polarity, so this should return False
        spectrum = Spectrum(mzml_fp, "test")
        spectrum.limit_polarity("negative")
        self.assertFalse(True in spectrum.to_use)

    def test_positive_polarity(self):
        # The entire spectrum is positive, so this should return False
        spectrum = Spectrum(mzml_fp, "test")
        spectrum.limit_polarity("positive")
        self.assertFalse(False in spectrum.to_use)

    def test_binning(self):
        spectrum = Spectrum(mzml_fp, "test")
        spectrum.limit_polarity("negative")
        spectrum.load_scans()

        # Create our own spectrum for testing purposes
        spectrum._masses = np.array([0.125, 0.1256, 0.130, 0.149, 0.1495])
        spectrum._intensities = np.array([50, 60, 100, 10, 50])

        # Apply binning
        spectrum.bin(0.005)

        # Hand calculated, like the sadistic person I am.
        self.assertTrue(np.array_equal(spectrum.masses, np.array([0.1253, 0.13, 0.14925])))
        self.assertTrue(np.array_equal(spectrum.intensities, np.array([55.0, 100, 30])))


    def test_infusion(self):
        spectrum = Spectrum(mzml_fp, "test")
        spectrum.limit_polarity("negative")
        spectrum.load_scans()

        spectrum._to_use = []

        class DuckyScan:
            def __init__(self, total_ion_count):
                self.TIC = total_ion_count

        tics = [10, 100, 200, 30, 40, 0, 0]

        for tic in tics:
            spectrum.scans.append(DuckyScan(tic))
            spectrum.to_use.append(True)

        spectrum._scans = np.array(spectrum.scans)
        spectrum._to_use = np.array(spectrum._to_use)

        spectrum.limit_infusion(1)

        self.assertTrue(np.array_equal(spectrum.to_use, np.array([False, True, True, False, False, False, False])))

    def test_reset(self):
        spectrum = Spectrum(mzml_fp, "test")
        spectrum.limit_polarity("negative")
        spectrum.load_scans()

        spectrum.reset()

        self.assertFalse(False in spectrum.to_use)

    def test_limit_spurious_peaks(self):
        # Dumb test, needs rewriting

        spectrum = Spectrum(mzml_fp, "test")
        spectrum.limit_polarity("positive")
        spectrum.load_scans()

        m_l = len(spectrum.masses)

        spectrum.remove_spurious_peaks()

        self.assertTrue(len(spectrum.masses) < m_l)

if __name__ == '__main__':
    unittest.main()