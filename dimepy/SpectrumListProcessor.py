# -*- coding: utf-8 -*-
# encoding: utf-8

import numpy as np
import warnings
import multiprocess
import scipy.stats as sc_stats
import bisect
import matplotlib.pyplot as plt
import pandas as pd
from sklearn.preprocessing import Imputer


class SpectrumListProcessor(object):
    def __init__(self, spectrum_list):
        '''

        :param spectrum_list:
        '''
        self.spectrum_list = spectrum_list
        self.mass_range = spectrum_list.get_mass_range()
        self._outlier_detected = False
        self._binned = False
        self._centered = False
        self._scaled = False
        self._value_imputated = False

    def to_list(self):
        return self.spectrum_list.to_list()

    def remove(self, spectrum):
        '''

        :param spectrum:
        :return:
        '''
        self.spectrum_list.remove(spectrum)

    def outlier_detection(self,
                          mad_threshold=3,
                          inplace=True,
                          plot=False,
                          results_path=None):
        '''

        :param mad_threshold:
        :param inplace:
        :param plot_path:
        :param results_path:
        :return:
        '''

        tics = [sum(s.intensities) for s in self.to_list()]

        mean_tic = np.nanmean(tics)
        mean_abs_dev = np.nanmean([abs(x - mean_tic) for x in tics])
        ad_f_m = [abs((x - mean_tic) / mean_abs_dev) for x in tics]

        outlier_spectrum = [
            s for i, s in enumerate(self.to_list())
            if ad_f_m[i] > mad_threshold
        ]

        if inplace is True:
            [self.remove(x) for x in outlier_spectrum]
            warnings.warn("Outlier detection removed: " +
                          ",".join([x.id for x in outlier_spectrum]))
            self._outlier_detected = True
        else:
            return outlier_spectrum

        if plot is True:
            plt.figure()
            plt.xlabel("Injection Order")
            plt.ylabel("Total Ion Count (TIC)")
            plt.scatter(
                [x._injection_order for x in self.to_list()],
                [sum(x.intensities) for x in self.to_list()],
                marker="o",
                color="b",
                label="Passed")
            plt.scatter(
                [x._injection_order for x in outlier_spectrum],
                [sum(x.intensities) for x in outlier_spectrum],
                marker="x",
                color="r",
                label="Outliers")
            plt.legend(loc="upper right", numpoints=1)
            plt.show()

    def binning(self, bin_size=0.25, statistic="mean", inplace=True, n_jobs=1):
        '''

        :param bin_size:
        :param statistic:
        :param inplace:
        :param n_jobs:
        :return:
        '''
        bins = np.arange(
            round(self.mass_range[0]),
            round(self.mass_range[1]),
            step=bin_size)

        def _make_bins():
            binned_masses = {b: [] for b in bins}

            for spectrum in self.to_list():
                for b in bins:
                    for m in spectrum.masses:
                        if m >= b and m < b + bin_size:
                            binned_masses[b].append(m)
            return binned_masses

        def _bin(spectrum):

            b_i, b_m, b_n = sc_stats.binned_statistic(
                spectrum.masses,
                spectrum.intensities,
                bins=bins,
                statistic=statistic)
            return b_m[:-1], b_i, spectrum.id

        bin_dicts = _make_bins()
        for idx, b in enumerate(bins):
            values = bin_dicts[b]
            if len(values) > 0:
                bins[idx] = sum(values) / len(values)

        pool = multiprocess.Pool(n_jobs)
        binned_spectra = pool.map_async(
            _bin, [spectrum for spectrum in self.to_list()]).get()
        pool.close()
        pool.join()

        if inplace is True:
            for result in binned_spectra:
                binned_masses, binned_intensities, id = result
                for spectrum in self.to_list():
                    if spectrum.id == id:
                        spectrum.masses = binned_masses
                        spectrum.intensities = binned_intensities
            self._binned = True

        else:
            warnings.warn("Non inplace binning yet to be implemented")

    def scale(self, method="mc", inplace=True):
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
        for spectrum in self.to_list():
            spectrum._normalise(method=method)

    def transform(self, method="nlog"):
        for spectrum in self.to_list():
            spectrum._transform(method=method)

    def value_imputation(self, method="basic", threshold=0.5, inplace=True):
        '''

        :param method:
        :param threshold:
        :param inplace:
        :return:
        '''

        def _remove_bins_by_threshold():
            df = self.spectrum_list.flatten_to_dataframe()

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

    def pandas_to_spectrum(self, df):
        '''

        :param df:
        :return:
        '''
        masses = df.columns
        for id, values in df.iterrows():
            intensities = values.values
            spectrum = [x for x in self.to_list() if x.id == id][0]
            spectrum.masses = masses
            spectrum.intensities = intensities

    def to_spectrumlist(self):
        '''

        :return:
        '''
        from SpectrumList import SpectrumList
        return SpectrumList(self.to_list())
