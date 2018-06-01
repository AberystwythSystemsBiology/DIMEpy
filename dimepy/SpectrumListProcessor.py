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
