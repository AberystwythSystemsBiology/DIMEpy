# Copyright (c) 2017-2019 Keiron O'Shea
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public
# License as published by the Free Software Foundation; either
# version 3 of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
#
# You should have received a copy of the GNU General Public
# License along with this program; if not, write to the
# Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
# Boston, MA 02110-1301 USA

import numpy as np
import math
from scipy.stats import binned_statistic
import pandas as pd

class SpectrumList:
    def __init__(self, spectrum_list: list = [], bin_width=0.01, statistic="mean"):
        self.bin_width = bin_width
        self.statistic = statistic
        if len(spectrum_list) > 1:
            self._spectrum_list = spectrum_list
            self.df = self.bin()
        else:
            # TODO: Throw exception if empty list passed
            exit(0)


    def bin(self):
        def _get_mass_range():
            mass_range = [x.mass_range for x in self._spectrum_list]
            min_mass = math.floor(np.min([x[0] for x in mass_range]))
            max_mass = math.ceil(np.max([x[1] for x in mass_range]) + self.bin_width)
            return min_mass, max_mass

        def _calculate_bins(bins):

            bin_dict = {x : [] for x in bins}

            intensities_dict = {}

            for spectrum in self._spectrum_list:
                masses = spectrum.masses
                intensities = spectrum.intensities

                binned_intensities, _, _ = binned_statistic(
                    masses,
                    intensities,
                    statistic=self.statistic,
                    bins=bins
                    )

                binned_masses, _, _ = binned_statistic(
                    masses,
                    masses,
                    statistic=self.statistic,
                    bins=bins
                )

                index = ~np.isnan(binned_intensities)


                binned_masses = binned_masses[index]

                for bin_indx, bin in enumerate(bins[:-1][index]):
                    bin_dict[bin].append(binned_masses[bin_indx])

                intensities_dict[spectrum.identifier] = binned_intensities

            return intensities_dict, bin_dict


        def calculate_masses(bin_dict):
            bins = []
            for bin in bin_dict:
                if bin_dict[bin] != []:
                    bins.append(np.mean(bin_dict[bin]))
                else:
                    bins.append(bin)
            return sorted(bins)

        min_mass, max_mass = _get_mass_range()
        bins = np.arange(min_mass, max_mass, step=self.bin_width)
        intensities_dict, bin_dict = _calculate_bins(bins)
        masses = calculate_masses(bin_dict)

        df = pd.DataFrame(intensities_dict).T
        df.columns = masses[:-1]
        df = df.loc[:, ~(pd.isna(df)).all(axis=0)]
        return df

    def transform(self, method: str ="log10"):

        def _trans(spec):
            if method.upper() == "LOG10":
                return np.log10(spec)
            elif method.upper() == "CUBE":
                return np.array(
                    [i**(1. / 3) for i in spec])
            elif method.upper() == "NLOG":
                return np.log(spec)
            elif method.upper() == "LOG2":
                return np.log2(spec)
            elif method.upper() == "GLOG":
                m = min(spec) / 10
                return np.log2(spec + np.sqrt(spec**2 + m**2)) / 2
            elif method.upper() == "SQRT":
                return np.array([sqrt(x) for x in spec])
            elif method.upper() == "IHS":
                return np.array([math.asinh(x) for x in spec])


        vals = self.df.values

        for i, x in enumerate(vals):
            vals[i] = _trans(x)

        self.df[:] = vals

    def normalise(self, method: str="tic"):
        def _normie(spec):
            if method.upper() == "TIC":
                sum_intensity = np.nansum(spec)
                normalised_intensities = np.array([(x / sum_intensity)
                                                   for x in spec]) * 1000
            elif method.upper() == "MEDIAN":
                median_intensity = np.nanmedian(spec)
                normalised_intensities = np.array(
                    [x - median_intensity for x in spec]) * 1000
            else:
                raise ValueError(
                    "%s is not a supported normalisation method" % method)
            return normalised_intensities

        vals = self.df.values

        for i, x in enumerate(vals):
            vals[i] = _normie(x)

        self.df[:] = vals

    def value_imputate(self, method: str="basic", threshold=0.5):

        def _remove_by_threshold():
            null_count = self.df.isnull().sum()
            _t = len(self.df.index.values) * threshold
            to_keep = null_count <= _t
            self.df[self.df.columns[to_keep.values]]


        def _apply_imputation():
            for identifier in self.df.index:
                i = self.df.ix[identifier].values
                if method.upper() == "BASIC":
                    filler = np.nanmin(i) / 2
                elif method.upper() == "MEAN":
                    filler = np.mean(i)
                elif method.upper() == "MIN":
                    filler = np.nanmin(i)
                elif method.upper() == "MEDIAN":
                    filler = np.nanmedian(i)
                else:
                    raise ValueError(
                        "%s is not a valid imputation method." % method
                        )
                i[np.isnan(i)] = filler
                self.df.ix[identifier] = i

        if method.upper() == "ALL":
            threshold = 0

        _remove_by_threshold()
        _apply_imputation()
        

    def tolist(self):
        return self._spectrum_list
