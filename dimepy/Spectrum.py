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
import warnings
import matplotlib.pyplot as plt
from operator import itemgetter
from Scans import Scans


class Spectrum(object):
    """A Spectrum class.

    This is a class that holds, manages, and manipulates mass spectrometry
    information.


    Parameters
    ----------

    fp : string, optional (default=None)
        Filepath to file of interest.

    id : string, optional (default=None)
        An identifier for the Spectrum.
        Note: if None is given one will be generated from the file path.

    polarity : string, optional (default=None)
        Polarity for fast polarity switching.

        - If "positive" then only positive ion scans will be read.
        - If "negative" then only negative ion scans will be read.
        - If None, all ion scans will be read.

    type : string, optional (default="peaks")
        Set what type of peaks you want to read into the Spectrum object.

        - If "peaks" then return all peaks of the spectrum.
        - If "centroidedPeaks" then return a centroided version of the
          profile spectrum.
        - If "reprofiledPeaks" then return a version a centroided spectrum.

    apex_mad : integer, optional (default=None)
        This method implements an 'automated' apex method akin to FIEmspro.

        - If None then return all scans as a single matrix.
        - If integer/float, then return all scans whose TIC are equal to or
          greater than the mean absolute deviation times the given value.

    snr_estimator : string, optional (default=None)
        Apply signal to noise filtering.

        - If None, then do not apply signal to noise filtering.
        - If "mad" then perform mean absolute deviation-based noise filtering.
        - If "mean" then perform mean-based noise filtering.
        - If "median" then perform median-based noise filtering.

    max_snr : float, optional (default=2.5)
        Value to use for signal to noise ratio filtering. If calculated snr is
        greater than or equal to the given value - then the signal will be
        removed.

    injection_order : integer, optional (default=None)
        The injection order of the Spectrum file.

    label : string, optional (default=None)
        Class label.

    Attributes
    ----------

    masses : array of shape = [n_masses] containing mass-to-ion data.
    intensities : array of shape = [n_intensities] containing intensity data.
    tic : integer calculated by summing the intensities.

    """

    def __init__(self,
                 fp=None,
                 id=None,
                 polarity=None,
                 type="peaks",
                 apex_mad=None,
                 snr_estimator=None,
                 max_snr=2.5,
                 injection_order=None,
                 label=None):

        self.fp = fp
        if polarity is not None:
            self.polarity = polarity.upper()
        else:
            self.polarity = polarity
        self.apex_mad = apex_mad
        self.max_snr = max_snr
        self.snr_estimator = snr_estimator
        self.type = type
        self.label = label

        self.masses = None
        self.intensities = None

        if self.fp is not None:
            self._load_from_file()
            if id is None:
                self._get_id_from_fp()
        if id is not None:
            self.id = id
        if injection_order is not None:
            self.injection_order = int(injection_order)
        else:
            self.injection_order = 0

        self._loaded = False

        self._normalised = False
        self._transformed = False
        self._scaled = False
        self._baseline_corrected = False

    @property
    def tic(self):
        return sum(self.intensities)

    def _get_id_from_fp(self):
        """Provides a spectrum identifier from the file_path attribute.

        If no id is passed, then the the filename (minus the file extension)
        will be used as the identifier.

        Notes
        -----

        Should not be called outside of this class.

        """
        self.id = os.path.splitext(os.path.basename(self.fp))[0]

    def baseline_correction(self, wsize=50, qtl=0.1, inplace=True):
        """Application of baseline correction over the spectrum intensities.

        Intensities are divided into equally spaced windows, where a local
        minimum intensity value (intervals are used to remove noise) is
        calculated as a baseline estimate for this region. The whole baseline
        is then computed through the use of linear interpolation via the
        central pairs of each interval.

        Parameters
        ----------

        wsize : integer (default=50)
            Window size to be used for indexing over intensities.

        qtl : float (default=0.1)
            A value to be used to calculate the lower quantile probabilitiy.

        inplace : boolean, optional (default=True)
            If False then return the corrected intensities array, else make the
            change within the object.

        Notes
        -----

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
            ymz = interv.iloc[::2, :].mean(axis=1).values

            ymz = interv.mean(axis=1).values
            ymz = np.insert(ymz, 0, 0)
            ymz = np.append(ymz, nmz)
            ymax = [y_max[0]
                    ] + y_max[0:] + [y_max[len(interv.index.values) - 1]]

            bl = np.interp(self.intensities, ymz, ymax)
            return bl

        if self._baseline_corrected is not True:
            bl = bc()
            baseline_corrected = self.intensities - bl
            if inplace is True:
                indx = baseline_corrected > 0
                self.masses = self.masses[indx]
                self.intensities = baseline_corrected[indx]
                self._baseline_corrected = True
            else:
                return baseline_corrected
        else:
            raise Exception("%s has already been baseline corrected" % self.id)

    def normalise(self, method="tic", inplace=True):
        """Application of normalisation over spectrum intensities.

        Normalisation aims to remove sources of variability within the spectrum

        Parameters
        ----------

        method : string, optional (default="tic")
            Method to use for normalisation.

            - If "tic" then apply TIC normalisation.
            - If "median" then apply normalisation by the median.

        inplace : boolean, optional (default=True)
            If False then return normalised intensities, else make the
            change within the object.

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
                raise ValueError(
                    "%s is not a supported normalisation method" % method)

        else:
            raise Exception("%s already normalised" % self.id)

        if inplace is True:
            self.intensities = normalised_intensities
            self._normalised = True
        else:
            return normalised_intensities

    def transform(self, method="log10", inplace=True):
        """Application of transformation over spectrum intensities.

        Transformation aims to make the data less skewed.

        Parameters
        ----------

        method : string, optional (default="log10")
            Method to use for transformation.

            - If "log10" then apply log10 transformation.
            - If "cube" then apply cube transformation.
            - If "nglog" then apply nlog transformation.
            - If "log2" then apply log2 transformation.
            - If "glog" then apply globalised log transformation.
            - If "sqrt" then apply square root-based transformation.
            - If "ihs" then apply inverse hyperbolic sine transformation.

        inplace : boolean, optional (default=True)
            If False then return normalised intensities, else make the
            change within the object.

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
                m = min(self.intensities) / 10
                transformed_intensities = np.log2(self.intensities + np.sqrt(
                    self.intensities**2 + m**2)) / 2
            elif method.upper() == "SQRT":
                transformed_intensities = np.array(
                    [sqrt(x) for x in self.intensities])
            elif method.upper() == "IHS":
                transformed_intensities = np.array(
                    [asinh(x) for x in self.intensities])

            else:
                raise ValueError(
                    "%s is not a supported transformation method" % method)
        else:
            raise Exception("%s already transformed" % self.id)

        if inplace is True:
            self.intensities = transformed_intensities
            self._transformed = True
        else:
            return transformed_intensities

    def _load_from_file(self):
        """Method for loading a Spectrum from file.

        Notes
        -----

        Should not be called outside of this class.

        """
        def __get_apex(scans):
            tics = scans.tics
            mad = np.mean(np.absolute(tics - np.mean(tics)))
            indx = tics >= mad * self.apex_mad
            scans.limiter(indx)

        scans = Scans(self.fp, self.snr_estimator, self.max_snr, self.type)

        if self.polarity is not None:
            indx = scans.polarities == self.polarity
            scans.limiter(indx)

        if self.apex_mad is not None:
            __get_apex(scans)

        masses, intensities = zip(*scans.scans)
        masses = np.concatenate(masses).ravel().tolist()
        intensities = np.concatenate(intensities).ravel().tolist()
        masses, intensities = zip(*sorted(zip(masses, intensities)))

        indx = np.array(intensities) != 0.0

        self.masses = np.array(masses)[indx]
        self.intensities = np.array(intensities)[indx]

    def plot(self, show=True, xlim=[], scaled=False, fp=None):
        """Method to visualise spectrum profile data using matplotlib.


        Parameters
        ----------

        show : boolean, optional (default=True)
            Whether or not to print the plot to screen.

        xlim : list, optional (default=[])
            Mass boundaries for plotting purposes, for example
            passing [50,500] would limit the plot to everything between those
            mass-to-ion values.

        scaled : boolean, optional (default=False)
             Whether or not to scale the intensities from 0 to 1.

        fp : string, optional (default=None)
            If provided, where to save the plot. Can accept all standard
            matplotlib outputs.

            - If None then do not save to file.

        """
        plt.figure()
        plt.title(self.id)
        plt.xlabel("Mass-to-ion (m/z)")
        plt.ylim(0, max(self.intensities))
        if self._normalised is True and self._transformed is False:
            plt.ylabel("Normalised Intensity")
        elif self._normalised is False and self._transformed is True:
            plt.ylabel("Transformed Intensity")
        elif self._normalised is True and self._transformed is True:
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
        if fp is not None:
            plt.savefig(fp)
        else:
            plt.show()
        plt.clf()
