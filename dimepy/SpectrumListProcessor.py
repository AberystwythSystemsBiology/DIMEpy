# -*- coding: utf-8 -*-

import numpy as np
import warnings
from pathos.multiprocessing import ProcessingPool as Pool
import scipy.stats as sc_stats
import bisect
import matplotlib.pyplot as plt
import pandas as pd
from sklearn.preprocessing import Imputer
from SpectrumList import SpectrumList
from copy import copy

class SpectrumListProcessor(SpectrumList):
    def __init__(self, spectrum_list):
        self._spectrum_list = copy(spectrum_list)
        self.mass_range = self._spectrum_list.get_mass_range()
        self._outlier_detected = False
        self._binned = False
        self._centered = False
        self._scaled = False
        self._value_imputated = False


    def outlier_detection(self,
                          threshold=3,
                          inplace=True,
                          plot=False):
        """Perform outlier detection.

        This method detects and removes sample outliers through the use of the
        mean absolute deviation over the total ion counts.

        Parameters
        ----------

        threshold : float, optional (default=3)

        plot : boolean, optional (default=False)

        inplace : boolean, optional (default=True)
            If False then return the corrected intensities array, else make the
            change within the object.

        """
        pass




    def binning(self,
                bin_size=0.25,
                int_statistic="median",
                mass_statistic="mean",
                inplace=True,
                n_jobs=1):
        """Perform mass-binning.

        Parameters
        ----------

        bin_size : float, optional (default=0.25)
            Mass-to-ion bin size to use for binning.

        int_statistic : string, optional (default="median")
            The method used to calculate the binned intensitiy value.

            - If mean, calculated as the mean all spectrum intensity values for
              a given bin.
            - If median, calculated as the median of all spectrum intensity values for
              a given bin.

        mass_statistic : string, optional (default="mean")
            The method used to calculate the binned mass value.

            - If mean, calculated as the mean all spectrum mass values for
              a given bin.
            - If median, calculated as the median of all spectrum mass values for
              a given bin.

        inplace : boolean, optional (default=True)
            If False then return the binned Spectrums, else make the
            change within the object.

        n_jobs : integer, optional (default=1)
            Number of threads to use to perform binning.

        """
        bins = np.arange(
            round(self.mass_range[0]),
            round(self.mass_range[1]),
            step=bin_size)

        def _make_bins():
            binned_masses = {b: [] for b in bins}
            for spectrum in self.tolist():
                sm = np.array(spectrum.masses)
                for b in bins:
                    m_bins = np.logical_and(sm >= b, sm <= (b + bin_size))
                    binned_masses[b].extend(sm[m_bins])
            return binned_masses

        def _bin(spectrum):
            b_i, b_m, b_n = sc_stats.binned_statistic(
                spectrum.masses,
                spectrum.intensities,
                bins=bins,
                statistic=int_statistic)
            return b_m[:-1], b_i, spectrum.id

        bin_dicts = _make_bins()

        for idx, b in enumerate(bins):
            values = bin_dicts[b]
            if len(values) > 0:
                if mass_statistic == "mean":
                    bins[idx] = np.mean(values).tolist()
                elif mass_statistic == "median":
                    bins[idx] = np.median(values).tolist()

        if len(self.tolist()) <= 2 or n_jobs == 1:
            binned_spectra = [_bin(s) for s in self.tolist()]
        else:
            pool = Pool(n_jobs)
            binned_spectra = pool.map(_bin, self.tolist())

        # TODO: This is awful. Needs to be reimplemented.
        if inplace is True:
            for result in binned_spectra:
                binned_masses, binned_intensities, id = result
                for spectrum in self.tolist():
                    if spectrum.id == id:
                        spectrum.masses = binned_masses
                        spectrum.intensities = binned_intensities
            self._binned = True
        else:
            masses = []
            ints = []
            ids = []
            for i in binned_spectra:
                binned_masses, binned_intensities, id = i
                masses = binned_masses
                ints.append(binned_intensities)
                ids.append(id)
            return pd.DataFrame(ints, columns=masses, index=ids)

    def scale(self, method="mc", inplace=True):
        """

        """
        def __mean_center(_i):
            return _i - np.mean(_i)

        def __pareto(_i):
            return __mean_center(_i) / np.sqrt(np.std(_i, dtype=np.float64))

        def __range(_i):
            return __mean_center(_i) / (np.max(_i) - np.min(_i))

        def __auto(_i):
            return __mean_center(_i) / np.std(_i)

        if self._scaled != True:
            df = self.to_spectrumlist().flatten_to_dataframe()
            for mass in df:
                intensities = df[mass].values
                if method.upper() == "MEAN":
                    scaled_intensities = __mean_center(intensities)
                elif method.upper() == "RANGE":
                    scaled_intensities = __range(intensities)
                elif method.upper() == "AUTO":
                    scaled_intensities = __auto(intensities)
                elif method.upper() == "PARETO":
                    scaled_intensities = __pareto(intensities)
                df[mass] = scaled_intensities

            self.pandas_to_spectrum(df)

    def normalise(self, method="tic"):
        for spectrum in self.tolist():
            spectrum._normalise(method=method)

    def transform(self, method="nlog"):
        for spectrum in self.tolist():
            spectrum._transform(method=method)

    def value_imputation(self, method="basic", threshold=0.5, inplace=True):
        def _remove_bins_by_threshold():
            df = self._spectrum_list.flatten_to_dataframe()

            nc = df.isnull().sum()
            _threshold = len(df.index.values) * threshold
            nc = nc[nc <= _threshold]
            return df[nc.index.values]

        def _value_imputation(df):
            if method.upper() == "KNN":
                imp = Imputer(axis=1)
                imputated = imp.fit_transform(df)
                df = pd.DataFrame(
                    imputated, columns=df.columns, index=df.index)
            else:
                for sample in df.index:
                    intensities = df.ix[sample].values
                    if method.upper() == "BASIC":
                        filler = np.nanmin(intensities) / 2
                    elif method.upper() == "MEAN":
                        filler = np.nanmean(intensities)
                    elif method.upper() == "MIN":
                        filler = np.nanmin(intensities)
                    elif method.upper() == "MEDIAN":
                        filler = np.nanmedian(intensities)
                    df.ix[sample] = df.ix[sample].replace(np.nan, filler)

            return df

        if method.upper() == "ALL":
            threshold = 0
            df = _remove_bins_by_threshold()
        else:
            df = _remove_bins_by_threshold()
            df = _value_imputation(df)

        if inplace is True:
            self.pandas_to_spectrum(df)
            self._value_imputated = True
        else:
            return df
