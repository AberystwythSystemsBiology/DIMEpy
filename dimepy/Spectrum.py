# -*- coding: utf-8 -*-

import os
from math import floor, asinh
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
from Scans import Scans

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
        label (str):

    Attributes:
        masses (np.array):
        intensities(np.array):
        _loaded (bool):
        _normalised (bool):
        _transformed (bool):
        _scaled (bool):
        _baseline_corrected (bool):
        _injection_order (int):
        masses (np.array):
        intensities (np.array):
        __raw_spectrum (pymzml.Reader):

    """

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
                 fp=None,
                 id=None,
                 polarity=None,
                 type="peaks",
                 apex=False,
                 apex_mad=2,
                 snr_estimator=None,
                 max_snr=2.5,
                 injection_order=None,
                 label=None):

        self.fp = fp
        if polarity != None:
            self.polarity = polarity.upper()
        else:
            self.polarity = polarity
        self.apex = apex
        self.apex_mad = apex_mad
        self.max_snr = max_snr
        self.snr_estimator = snr_estimator
        self.type = type
        self.label = label

        if self.fp != None:
            self._load_from_file()
            if id is None:
                self._get_id_from_fp()
        if id != None and self.id != None:
            self.id = id

        if injection_order is not None:
            self._injection_order = int(injection_order)
        else:
            self._injection_order = 0



    def _get_id_from_fp(self):
        """Provides a spectrum identifier from the file_path attribute.

        If no id is passed, then the the filename (minus the file extension)
        will be used as the identifier.

        Note:
            Should not be ran outside of this class.
        """
        self.id = os.path.splitext(os.path.basename(self.fp))[0]


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
        if self._baseline_corrected != True:
            bl = bc()
            baseline_corrected = self.intensities - bl
            if inplace == True:
                indx = baseline_corrected > 0
                self.masses = self.masses[indx]
                self.intensities = baseline_corrected[indx]
            else:
                return baseline_corrected
        else:
            raise Exception("%s has already been baseline corrected" % self.id)

    def normalise(self, method="tic", inplace=True):
        """Application of normalisation over the spectrum intensities.

        Normalisation aims to remove sources of variability within the spectrum.

        Args:
            method (str):
            inplace (bool): If False then return the corrected baseline array,
                else make the change within the object.
        """

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
                raise ValueError("%s is not a supported normalisation method" % method)

        else:
            raise Exception("%s already normalised" % self.id)

        if inplace == True:
            self.intensities = normalised_intensities
            self._normalised = True
        else:
            return normalised_intensities

    def transform(self, method="log10", inplace=True):
        """

        Args:
            method (str):
            inplace (bool):

        """

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
                min = min(self.intensities) / 10
                transformed_intensities = np.log2(self.intensities + np.sqrt(self.intensities**2 + min**2)) / 2
            elif method.upper() == "SQRT":
                transformed_intensities = np.array([sqrt(x) for x in self.intensities])
            elif method.upper() == "IHS":
                transformed_intensities = np.array([asinh(x) for x in self.intensities])
            else:
                raise ValueError("%s is not a supported transformation method" % method)
        else:
            raise Exception("%s already transformed" % self.id)

        if inplace is True:
            self.intensities = transformed_intensities
            self._transformed = True
        else:
            return transformed_intensities

    def _load_from_file(self):
        def __get_apex(scans):
            tics = scans.tics
            mad = np.mean(np.absolute(tics - np.mean(tics)))
            indx = tics >= mad*self.apex_mad
            scans.limiter(indx)

        scans = Scans(self.fp, self.snr_estimator, self.max_snr, self.type)

        if self.polarity != None:
            indx = scans.polarities == self.polarity
            scans.limiter(indx)

        if self.apex == True:
            __get_apex(scans)

        masses, intensities = zip(*scans.scans)
        masses = np.concatenate(masses).ravel().tolist()
        intensities = np.concatenate(intensities).ravel().tolist()

        masses, intensities = zip(*sorted(zip(masses, intensities)))
        self.masses = np.array(masses)
        self.intensities = np.array(intensities)

    def plot(self, show=True, xlim=[], scaled=False, file_path=None):
        """Method to visualise spectrum profile data using matplotlib.

        Args:
            show (boolean): Whether or not to print the plot to screen.
            xlim (list): Mass boundaries for plotting purposes, for example
                passing [50,500] would limit the plot to everything between those
                values.
            scaled (boolean): Whether or not to scale the intensities from 0 to 1.
            file_path (str): If provided, where to save the plot. Can accept all
                standard matplotlib outputs.

        """
        plt.figure()
        plt.title(self.id)
        plt.xlabel("Mass-to-ion (m/z)")
        plt.ylim(0, max(self.intensities))
        if self._normalised == True and self._transformed == False:
            plt.ylabel("Normalised Intensity")
        elif self._normalised == False and self._transformed == True:
            plt.ylabel("Transformed Intensity")
        elif self._normalised == True and self._transformed == True:
            plt.ylabel("Processed Intensity")
        else:
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
