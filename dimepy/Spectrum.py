# -*- coding: utf-8 -*-

import os
from math import floor
from collections import defaultdict
from copy import copy
import numpy as np
import pandas as pd
import operator
import pymzml
from scipy.stats.mstats import mquantiles
from scipy.sparse import csc_matrix, eye, diags
from scipy.sparse.linalg import spsolve
from scipy.interpolate import interp1d
import warnings
import matplotlib.pyplot as plt
from operator import itemgetter

np.seterr("ignore")

class Spectrum(object):
    """This is a class that holds, manages, and manipulates mass spectrometry
    data/information.

    Args:
        file_path (str):
        id (str):
        polarity (str):
        type (str):
        min_mz (float):
        max_mz (float):
        apex (bool):
        apex_mad (int):
        snr_estimator (str):
        injection_order (int):

    Attributes:
        _loaded (bool):
        _normalised (bool):
        _scaled (bool):
        _baseline_corrected (bool):
        _injection_order (int):
        masses (np.array):
        intensities (np.array):
        __raw_spectrum (pymzml.Reader):

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
                 type="peaks",
                 min_mz=50.0,
                 max_mz=5000,
                 apex=False,
                 apex_mad=2,
                 snr_estimator="median",
                 max_snr=2.5,
                 injection_order=None):

        self.file_path = file_path
        if id is None and self.file_path != None:
            self._get_id_from_fp()
        else:
            self.id = id

        if injection_order is not None:
            self._injection_order = int(injection_order)
        else:
            self._injection_order = 0
        self.apex = apex
        self.apex_mad = apex_mad
        self.max_snr = max_snr
        self.snr_estimator = snr_estimator
        self.type = type
        self.min_mz = min_mz
        self.max_mz = max_mz
        self.polarity = polarity
        if self.file_path != None:
            self._load_from_file()
        self._normalised = False
        self._transformed = False
        self._baseline_corrected = False
        self._scaled = False

    def _get_id_from_fp(self):
        """Provides a spectrum identifier from the file_path attribute.

        If no id is passed, then the the filename (minus the file extension)
        will be used as the identifier.
        """
        self.id = os.path.splitext(os.path.basename(self.file_path))[0]

    def baseline_correction(self, wsize=50, qtl=0.1, inplace=True):
        """Application of baseline correction over the spectrum intensities.

        Intensities are divided into equally spaced windows, where a local minimum
        intensity value (intervals are used to remove noise) is calculated as a baseline
        estimate for this region. The whole baseline is then computed through the use
        of linear interpolation via the central pairs of each interval.

        Args:
            wsize (int): Window size to be used for indexing over intensities.
            qtl (float): A value to be used to calculate the lower quantile probabilitiy.
            inplace (bool): If False then return the corrected baseline array,
                else make the change within the object.

        Note:
            Python implementation of onebc.R from FIEmspro.
        """

        def bc():
            # TODO: Pythonise this.
            x = self.intensities
            nmz = len(x)
            addini = floor(nmz % wsize / 2)
            i = np.arange(addini + wsize, nmz - wsize, wsize)
            i = np.insert(i, 0, 0)
            j = np.arange(addini + wsize, nmz - wsize, wsize)
            j = np.append(j, nmz)
            interv = pd.DataFrame(zip(*[i, j]))
            interv = interv.append(
                interv.iloc[1:] - floor(wsize / 2), ignore_index=True)
            interv = interv.sort_values(0, ascending=True)
            y_max = []
            for _, k in interv.iterrows():
                k_b = x[int(k[0]):int(k[1])]
                y_max.append(mquantiles(k_b, prob=qtl)[0])
            #ymz = interv.iloc[::2, :].mean(axis=1).values

            ymz = interv.mean(axis=1).values
            ymz = np.insert(ymz, 0, 0)
            ymz = np.append(ymz, nmz)
            ymax = [y_max[0]
                    ] + y_max[0:] + [y_max[len(interv.index.values) - 1]]

            bl = np.interp(self.intensities, ymz, ymax)
            return bl

        bl = bc()
        baseline_corrected = self.intensities - bl
        if inplace == True:
            indx = baseline_corrected > 0
            self.masses = self.masses[indx]
            self.intensities = baseline_corrected[indx]
        else:
            return baseline_corrected

    def _normalise(self, method="tic", inplace=True):
        """Application of normalisation over the spectrum intensities.

        Args:
            method (str):
            inplace (bool): If False then return the corrected baseline array,
                else make the change within the object.
        """"

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

        if inplace == True:
            self.intensities = normalised_intensities
            self._normalised = True
        else:
            return normalised_intensities

    def _transform(self, method="log10", inplace=True):
        """"

        Args:
            method (str):
            inplace (bool):

        """"

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
        def __from_mzml():
            def __gen_reader():
                return pymzml.run.Reader(
                    self.file_path,
                    extraAccessions=[["MS:1000129", ["value"]],
                                     ["MS:1000130", ["value"]]])

            def __get_scans_of_interest():
                reader = __gen_reader()
                scans = []
                for scan_number, scan in enumerate(reader):
                    if self.polarity != None:
                        if scan.get(self.polarity_dict[self.polarity.upper()]
                                    ) != None:
                            scans.append(scan_number)
                    else:
                        scans.append(scan_number)

                reader = __gen_reader()
                if self.apex == True:
                    '''
                        The following method is taken from FIEmspro, with a small
                        amount of modification for automated selection of the
                        infusion profile scan-ranges.
                    '''
                    tics = []
                    for scan_number, scan in enumerate(reader):
                        if scan_number in scans:
                            tic = sum(zip(*scan.peaks)[1])
                            tics.append([scan_number, tic])
                    mad = np.mean(
                        np.absolute(zip(*tics)[1] - np.mean(zip(*tics)[1])))

                    peak_scans = [
                        x for x in tics if x[1] >= (self.apex_mad * mad)
                    ]
                    # I've noticed that some profiles have a strange overflow at the end
                    # of the run, so I will look for those and discard here.
                    peak_range_mad = np.mean(
                        np.absolute(
                            zip(*peak_scans)[0] - np.mean(zip(*peak_scans)[0]))
                    )
                    scans = [
                        x[0] for x in peak_scans if x[0] >= peak_range_mad
                    ]

                return scans

            def __get_scan(scan_range):
                masses = []
                intensities = []
                reader = __gen_reader()
                for scan_number, scan in enumerate(reader):
                    if scan["ms level"] != None and scan_number in scan_range:
                        m, ints = zip(*scan.peaks)
                        m, ints = np.array(m), np.array(ints)
                        indx = np.logical_and(
                            np.array(m) >= self.min_mz,
                            np.array(m) <= self.max_mz)
                        m = m[indx]
                        ints = ints[indx]
                        if len(m) > 0:
                            if self.snr_estimator != None:
                                sn_r = np.divide(
                                    ints,
                                    scan.estimatedNoiseLevel(
                                        mode=self.snr_estimator))
                                gq = sn_r >= self.max_snr
                                m = m[gq]
                                ints = ints[gq]
                            masses.extend(m)
                            intensities.extend(ints)

                masses, intensities = zip(*sorted(zip(masses, intensities)))
                masses, intensities = np.array(masses), np.array(intensities)
                non_z = [intensities != 0]
                return masses[non_z], intensities[non_z]

            scan_range = __get_scans_of_interest()
            self.masses, self.intensities = __get_scan(scan_range)

        if self.file_path.upper().endswith("MZML"):
            __from_mzml()
        elif self.file_path.upper().endswith("CSV"):
            print "CSV not implemented"
        elif self.file_path.upper().endswith("PKL"):
            print "PKL not implemented"
        else:
            raise Exception()

    def plot(self, show=True, xlim=[], scaled=False, file_path=None):
        """

        Args:
            show (boolean):
            xlim (list):
            scaled (boolean):
            file_path (str):

        """
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
