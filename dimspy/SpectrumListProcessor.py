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

    def outlier_detection(self, mad_threshold=3, inplace=True, plot=False, results_path=None):
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

        outlier_spectrum = [s for i, s in enumerate(
            self.to_list()) if ad_f_m[i] > mad_threshold]

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
            plt.scatter([x._injection_order for x in self.to_list()],
                        [sum(x.intensities) for x in self.to_list()],
                        marker="o", color="b", label="Passed")
            plt.scatter([x._injection_order for x in outlier_spectrum],
                        [sum(x.intensities) for x in outlier_spectrum],
                        marker="x", color="r", label="Outliers")
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
        def _bin(spectrum):
            bins = np.arange(
                self.mass_range[0], self.mass_range[1], step=bin_size)

            b_i, b_m, b_n = sc_stats.binned_statistic(spectrum.masses,
                                                      spectrum.intensities,
                                                      bins=bins,
                                                      statistic=statistic)

            return [bins[:-1], b_i, spectrum.id]

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

    def scale(self, method="MC", inplace=True, n_jobs=1):
        def _mean_center(spectrum):
            mean_intensity = np.nanmean(spectrum.intensities)
            return np.array([x - mean_intensity for x in spectrum.intensities])

        def _scaler(data):
            spectrum, method = data
            mean_centered = _mean_center(spectrum)
            if method.upper() == "AUTO":
                variance = np.var(mean_centered)
                scaled_intensities = np.array(
                    [x / variance for x in mean_centered])
            elif method.upper() == "RANGE":
                r = np.nanmax(spectrum.intensities) - \
                    np.nanmin(spectrum.intensities)
                scaled_intensities = np.array([x / r for x in mean_centered])
            elif method.upper() == "PARETO":
                pareto = np.std(spectrum.intensities)**(1 / 2)
                scaled_intensities = np.array(
                    [x / pareto for x in mean_centered])
            elif method.upper() == "MC":
                scaled_intensities = mean_centered
            return scaled_intensities, spectrum.id

        pool = multiprocess.Pool(n_jobs)

        scaled_spectrum = pool.map_async(_scaler, [[spectrum, method] for spectrum in
                                                   self.to_list()])
        scaled_spectrum = scaled_spectrum.get()

        pool.close()
        pool.join()

        for result in scaled_spectrum:
            scaled_intensities, id = result
            if inplace is True:
                for spectrum in self.to_list():
                    if spectrum.id == id:
                        spectrum.intensities = scaled_intensities
                self._scaled = True
            else:
                warnings.warn("Non-inplace centering not implemented!")

    def normalise(self, method="tic"):
        for spectrum in self.to_list():
            spectrum._normalise(method=method)

    def transform(self, method="nlog"):
        for spectrum in self.to_list():
            spectrum._normalise(method=method)

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
            df = df[nc[nc < (len(df.index.values) * threshold)].index.values]
            return df

        def _value_imputation(df):
            if method.upper() == "KNN":
                imp = Imputer(axis=0)
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
            df = self.spectrum_list.flatten_to_dataframe()
            df = df.loc[:, df.isnull().mean() < 1]
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
