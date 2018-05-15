# -*- coding: utf-8 -*-

import os
import numpy as np
import operator
import pymzml
from scipy.sparse import csc_matrix, eye, diags
from scipy.sparse.linalg import spsolve
import warnings
import matplotlib.pyplot as plt


class Spectrum(object):
    """Summary of class here.

    Longer class information...

    Attributes:
        something: description of said thing.

    """
    polarity_dict = {"POSITIVE": "MS:1000130", "NEGATIVE": "MS:1000129"}

    _loaded = False

    _normalised = False
    _transformed = False
    _scaled = False
    _baseline_corrected = False

    _injection_order = None

    masses = np.array([])
    intensities = np.array([])

    __raw_spectrum = None

    def __init__(self,
                 file_path=None,
                 id=None,
                 polarity=None,
                 parameters=None,
                 injection_order=None):
        """Inits a Spectrum.

        Forms a Spectrum object, contained masses and intensities..

        Args:
            file_path:
            id:
            polarity:
            parameters:
            injection_order:

        """
        self.file_path = file_path
        if id is None:
            self._get_id_from_fp()
        else:
            self.id = id

        if parameters is None:
            self.parameters = {
                "MS1 Precision": 1e-3,
                "MSn Precision": 1e-3,
                "Measured Precision": 1e-3,
                "Scan Range": "apex",
                "Peak Type": "centroided"
            }
        else:
            self.parameters = parameters

        if injection_order is not None:
            self._injection_order = int(injection_order)
        else:
            warnings.warn(
                "Injection order information is required if you want to use outlier detection!"
            )

        self.polarity = polarity
        if self.file_path != None:
            self._load_from_file()
        self._normalised = False
        self._transformed = False
        self._baseline_corrected = False
        self._scaled = False

    def _get_id_from_fp(self):
        """Description of method.

        Longer description of method.
        """
        self.id = os.path.splitext(os.path.basename(self.file_path))[0]

    # Over-egging the WhittakerSmoothing, need to take a look.

    def baseline_correction(self, inplace=True, max_iterations=2, lambda_=100):
        """Description of method.

        Longer description of method.
        """
        warnings.warn("This is currently in development...")

        def _WhittakerSmooth(intensities_copy, ones):
            intensities_copy = np.matrix(intensities_copy)
            diag_eye = eye(intensities_copy.size, format="csc")
            diag = diag_eye[1:] - diag_eye[:-1]
            sparse = diags(
                ones, 0, shape=(intensities_copy.size, intensities_copy.size))

            csc_A = csc_matrix(sparse + (lambda_ * diag.T * diag))
            csc_B = csc_matrix(sparse * intensities_copy.T)
            return np.array(spsolve(csc_A, csc_B))

        def _AirPLS():
            intensities_copy = self.intensities
            ones = np.ones(self.intensities.shape[0])
            for index in range(0, max_iterations):
                whittaker_smoothed = _WhittakerSmooth(intensities_copy, ones)
                smoothed_intensities = (intensities_copy - whittaker_smoothed)
                smoothed_sum = np.abs(
                    smoothed_intensities[smoothed_intensities < 0].sum())
                if (smoothed_sum < 0.001 * (abs(intensities_copy)).sum()
                        or index == max_iterations):
                    break
                ones[smoothed_intensities >= 0] = [0]
                ones[smoothed_intensities < 0] = np.exp(
                    index * (smoothed_intensities[smoothed_intensities < 0]) /
                    smoothed_sum)
                ones[0] = np.exp(
                    index *
                    (smoothed_intensities[smoothed_intensities < 0]).max() /
                    smoothed_sum)
                ones[-1] = ones[0]
            return smoothed_intensities

        if self._baseline_corrected is True:
            warnings.warn(
                "It seems like this spectrum has already been baseline corrected!"
            )

        calculated_baseline = _AirPLS()

        baseline_corrected_intensities = []
        baseline_corrected_masses = []

        for index, intensity in enumerate(
                self.intensities - calculated_baseline):
            if intensity > 0:
                baseline_corrected_intensities.append(intensity)
                baseline_corrected_masses.append(self.masses[index])

        baseline_corrected_intensities = np.array(
            baseline_corrected4_intensities)
        baseline_corrected_masses = np.array(baseline_corrected_masses)

        if inplace is True:
            self.intensities = baseline_corrected_intensities
            self.masses = baseline_corrected_masses
            self._baseline_corrected = True
        else:
            return baseline_corrected_masses, baseline_corrected_intensities

    def _normalise(self, method="tic"):
        '''

        :param method:
        :param inplace:
        :return:
        '''

        if self._normalised is False:
            if method.upper() == "TIC":
                sum_intensity = np.nansum(self.intensities)
                normalised_intensities = np.array([(
                    x / sum_intensity) for x in self.intensities]) * 1000
            elif method.upper() == "MEDIAN":
                median_intensity = np.nanmedian(self.intensities)
                normalised_intensities = np.array(
                    [x - median_intensity for x in self.intensities]) * 1000
            else:
                normalised_intensities = self.intensities
        else:
            warnings.warn("Warning: %s already normalised, ignoring" % self.id)
            normalised_intensities = self.intensities

        self.intensities = normalised_intensities
        self._normalised = True

    def _transform(self, method="log10", inplace=True):
        '''

        :param method:
        :param inplace:
        :return:
        '''

        if self._transformed is False:
            if method.upper() == "LOG10":
                transformed_intensities = np.log10(self.intensities)
            elif method.upper() == "CUBE":
                transformed_intensities = np.array(
                    [i**(1. / 3) for i in self.intensities])
            elif method.upper() == "NLOG":
                transformed_intensities = np.log(self.intensities)
            elif method.upper() == "LOG2":
                transformed_intensities = np.log2(self.intensities)
            elif method.upper() == "GLOG":

                def __lognorm(x, min):
                    return np.log2((x + np.sqrt(x**2 + min**2)) / 2)

                transformed_intensities = __lognorm(self.intensities,
                                                    min(self.intensities) / 10)
            else:
                pass
        else:
            warnings.warn("Warning: %s already normalised, ignoring" % self.id)
            transformed_intensities = self.intensities

        if inplace is True:
            self.intensities = transformed_intensities
            self._transformed = True
        else:
            return transformed_intensities

    def _load_from_file(self):
        if self.file_path.upper().endswith("MZML"):
            self._from_mzml()
        elif self.file_path.upper().endswith("CSV"):
            print "CSV not implemented"
        elif self.file_path.upper().endswith("PKL"):
            print "PKL not implemented"
        else:
            raise Exception()

    def _from_mzml(self):
        def __polarity():
            reader = pymzml.run.Reader(
                self.file_path,
                extraAccessions=[["MS:1000129", ["value"]],
                                 ["MS:1000130", ["value"]]])

            polarity_scans = []
            for scan_number, scan in enumerate(reader):
                if scan.get(self.polarity_dict[self.polarity.upper()]) != None:
                    polarity_scans.append(scan_number)
            return polarity_scans

        def __return_apex(polarity_scans):
            reader = pymzml.run.Reader(
                self.file_path,
                MSn_Precision=self.parameters["MSn Precision"],
                MS1_Precision=self.parameters["MS1 Precision"])
            tic_scans = []
            for scan_number, scan in enumerate(reader):
                if scan_number in polarity_scans:
                    tic_scans.append([scan["total ion current"], scan_number])
            tics = [x[0] for x in tic_scans]
            mad = np.mean(np.absolute(tics - np.mean(tics))) * 3

            __get_spectrum([x[1] for x in tic_scans if x[0] > mad])

        def __get_spectrum(scans):
            reader = pymzml.run.Reader(
                self.file_path,
                MS1_Precision=self.parameters["MS1 Precision"],
                MSn_Precision=self.parameters["MSn Precision"])

            self.__raw_spectrum = pymzml.spec.Spectrum(
                measuredPrecision=self.parameters["Measured Precision"])

            for scan_number, scan in enumerate(reader):
                if scan_number in scans:
                    self.__raw_spectrum += scan

            if self.parameters["Peak Type"] == "centroided":
                spectrum = [[m, i]
                            for m, i in self.__raw_spectrum.centroidedPeaks]
            elif self.parameters["Peak Type"] == "reprofiled":
                spectrum = [[m, i]
                            for m, i in self.__raw_spectrum.reprofiledPeaks]
            else:
                spectrum = [[m, i] for m, i in self.__raw_spectrum.peaks]

            spectrum = sorted(spectrum, key=operator.itemgetter(0))

            self.masses = np.array([x[0] for x in spectrum])
            self.intensities = np.array([x[1] for x in spectrum])

            self._loaded = True

        polarity_scans = __polarity()

        if self.parameters["Scan Range"] == "apex":
            __return_apex(polarity_scans)
        else:
            __get_spectrum(polarity_scans)

    def plot(self, show=True, xlim=[], scaled=False, file_path=None):
        plt.figure()
        plt.title(self.id)
        plt.xlabel("Mass-to-ion (m/z)")
        plt.ylim(0, max(self.intensities))
        plt.ylabel("Intensity")

        if xlim == []:
            xlim = [min(self.masses), max(self.masses)]
            plt.ylim(0, max(self.intensities))
        plt.xlim(xlim)
        if scaled is False:
            plt.ylim(0, max(self.intensities))
            plt.plot(self.masses, self.intensities)
            plt.ylabel("Intensity")
        else:
            scaled_intensities = self.normalise(method="tic", inplace=False)
            plt.plot(self.masses, scaled_intensities)
            plt.ylim(0, max(scaled_intensities))
            plt.ylabel("Scaled Intensity")
        plt.tight_layout()
        if file_path is not None:
            plt.savefig(file_path)
        else:
            plt.show()
        plt.clf()
