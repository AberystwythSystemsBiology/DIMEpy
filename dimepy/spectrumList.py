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
from scipy.stats import binned_statistic
from .spectrum import Spectrum
import math
from typing import Tuple, List
from .utils import bin_masses_and_intensities
import csv


class SpectrumList:

    def __init__(self):
        self._list = []

        self.binned = False
        self.value_imputated = False
        self.normalised = False
        self.transformed = False

        self._global_masses = False

    def append(self, spectrum: Spectrum):
        """
        Method to append a Spectrum to the SpectrumList.

        Arguments:
            spectrum (Spectrum): A Spectrum object.
        """
        if type(spectrum) == Spectrum:
            self._list.append(spectrum)
        else:
            raise ValueError("SpectrumList only accepts Spectrum objects.")

    def detect_outliers(self, threshold: float = 1, verbose: bool = False):
        """
        Method to locate and remove outlier spectrum using MAD.

        Arguments:
            threshold (int): Threshold for MAD outlier detection.
        """

        def _get_tics() -> Tuple[np.array, np.array]:
            tics = []

            for spec in self._list:
                tics.append(np.sum(spec.intensities))

            return np.array(tics)

        def _calculate_mad(tics: np.array) -> float:
            return np.median(np.abs(tics - np.median(tics)))

        def _get_mask(tics: np.array, mad: float) -> np.array:
            tics = tics[:, None]
            median = np.median(tics, axis=0)
            diff = np.sum((tics - median)**2, axis=-1)
            diff = np.sqrt(diff)

            med_abs_deviation = np.median(diff)

            modified_z_score = 0.6745 * diff / med_abs_deviation

            return modified_z_score <= threshold

        tics = _get_tics()
        mad = _calculate_mad(tics)
        to_keep = _get_mask(tics, mad)

        _list = np.array(self._list)

        if verbose:
            print("Detected Outliers: %s" %
                  ";".join([x.identifier for x in _list[~to_keep]]))

        self._list = _list[to_keep].tolist()

    def bin(self, bin_width: float = 0.5, statistic: str = "mean"):
        """
        Method to conduct mass binning to nominal mass and mass spectrum
        generation across a SpectrumList.

        Arguments:
            bin_width (float): The mass-to-ion bin-widths to use for binning.
            statistic (str): The statistic to use to calculate bin values.
        """

        def _get_global_mass_range() -> Tuple[float, float]:
            mass_ranges = [s.mass_range for s in self._list]

            min_mass = min([_min for _min, _max in mass_ranges]) - bin_width
            max_mass = max([_max for _min, _max in mass_ranges]) + bin_width

            return min_mass, max_mass

        def _get_global_bins(min_mass: float, max_mass: float):
            bins = np.arange(min_mass, max_mass, step=bin_width)

            bin_dict = {x: [] for x in bins}
            intensities = []

            for spec in self._list:
                m = spec.masses
                i = spec.intensities

                binned_i, _, _ = binned_statistic(m,
                                                  i,
                                                  statistic=statistic,
                                                  bins=bins)

                binned_m, _, _ = binned_statistic(m,
                                                  m,
                                                  statistic=statistic,
                                                  bins=bins)

                index = ~np.isnan(binned_i)

                binned_m = binned_m[index]
                binned_i = binned_i[index]

                for indx, b in enumerate(bins[np.where(index == True)]):
                    bin_dict[b].append(binned_m[indx])

                intensities.append([binned_i, index])

            return bin_dict, intensities

        def _get_masses(bin_dict) -> np.array:
            bins = []
            for b, bi in bin_dict.items():
                if bi != []:
                    bins.append(np.mean(bi))
                else:
                    bins.append(b)
            return np.array(bins)[:-1]

        min_mass, max_mass = _get_global_mass_range()
        bin_dict, intensities = _get_global_bins(min_mass, max_mass)
        self._global_masses = _get_masses(bin_dict)

        # Apply to spectrum objects
        for index, (intensities, mass_index) in enumerate(intensities):
            s = self._list[index]
            s._masses = self._global_masses[mass_index]
            s._intensities = intensities

        self.binned = True

    def value_imputate(self, method: str = "min",
                       threshold: float = 0.5) -> None:
        """
        A method to deploy value imputation to our spectrum list.

        Arguments:
            method (str): Method to use for value imputation.
            threshold (float): Number of samples an intensity needs to be
            present in to be taken forward for imputation.
        """

        def _extend_spectrum():
            for spec in self._list:
                m = spec.masses
                i = spec.intensities

                is_in = np.intersect1d(self._global_masses,
                                       m,
                                       return_indices=True)[1]

                # Empty intensities
                exp_i = np.empty(self._global_masses.shape)
                # Replace with np.nans
                exp_i[:] = np.nan
                # Put the mass masked data in.
                exp_i[is_in] = i

                spec._masses = self._global_masses
                spec._intensities = exp_i

        def _determine_to_keep():
            global_intensities = np.array([s.intensities for s in self._list])
            global_intensities[np.isnan(global_intensities)] = 0.0
            non_zeros = np.count_nonzero(global_intensities, axis=0)
            return non_zeros >= global_intensities.shape[0] * threshold

        def _apply_to_keep(to_keep):
            self._global_masses[to_keep]

            for spec in self._list:
                spec._masses = spec.masses[to_keep]
                spec._intensities = spec.intensities[to_keep]

        def _apply_imputation():
            for s in self._list:
                i = s.intensities

                if method.upper() == "BASIC":
                    filler = np.nanmin(i) / 2
                elif method.upper() == "MEAN":
                    filler = np.mean(i)
                elif method.upper() == "MIN":
                    filler = np.nanmin(i)
                elif method.upper() == "MEDIAN":
                    filler = np.nanmedian(i)
                else:
                    raise ValueError("%s is not a valid imputation method." %
                                     method)

                i[np.isnan(i)] = filler

                s._intensities = i

        if self.binned:
            _extend_spectrum()
            to_keep = _determine_to_keep()
            _apply_to_keep(to_keep)
            _apply_imputation()
            self.value_imputated = True

        else:
            raise ValueError(
                "This only works where the SpectrumList has been binned.")

    def normalise(self, method: str = "tic") -> None:
        """
        Method to conduct sample independent intensity normalisation.

        Arguments:
            method (str): The normalisation method to use.
        """

        def _normie(spec: Spectrum):
            i = spec.intensities

            if method.upper() == "TIC":
                spec._intensities = np.divide(i, np.sum(i)) * 1000
            elif method.upper() == "MEDIAN":
                spec._intensities = i - np.median(i) * 1000
            elif method.upper() == "MEAN":
                spec._intensities = i - np.mean(i) * 1000
            else:
                raise ValueError("%s is not a valid normalisation method" %
                                 (method))

        if self.normalised:
            for spec in self._list:
                _normie(spec)

            self.normalised = True
        else:
            raise ValueError(
                "It looks like you've already normalised this data.")

    def transform(self, method: str = "log10") -> None:
        """
        Method to conduct sample independent intensity transformation.

        Arguments:
            method (str): The transformation method to use.
        """

        def _transform(spec: Spectrum):
            i = spec.intensities

            if method.upper() == "LOG10":
                spec._intensities = np.log10(i)
            elif method.upper() == "CUBE":
                spec._intensities = i**(1. / 3)
            elif method.upper() == "NLOG":
                spec._intensities = np.log(i)
            elif method.upper() == "LOG2":
                spec._intensities = np.log2(i)
            elif method.upper() == "GLOG":
                m = np.min(i) / 10
                spec._intensities = np.log2(i + np.sqrt(i**2 + m**2)) / 2
            elif method.upper() == "SQRT":
                spec._intensities = np.sqrt(i)
            elif method.upper() == "IHS":
                spec._intensities = np.array([math.asinh(x) for x in i])

        if self.transformed:
            for spec in self._list:
                _transform(spec)
            self.transformed = True
        else:
            raise ValueError(
                "It looks like you've already transformed this data.")

    def to_csv(self,
               fp: str,
               sep: str = ",",
               output_type: str = "metaboanalyst"):
        """
        Method to export the spectrum list.

        Arguments:
            fp (str): Filepath to export the file to.
            sep (str): Seperator to use for file export
            output_type (str): Output type.
        """

        def _to_metaboanalyst():
            _output = []

            for s in self._list:
                _samp = np.array([s.identifier, s.stratification])
                _samp = np.append(_samp, np.array([s.masses, s.intensities]).T)
                _samp = _samp.reshape((s.intensities.shape[0] + 1, 2))

                _output.append(_samp)
            _output = np.array(_output)

            _output = np.hstack(_output)
            np.savetxt(fp, _output, delimiter=sep, fmt="%s")

        def _to_matrix():
            _output = np.ndarray(
                (len(self._list[0].masses) + 1, len(self._list) + 1),
                dtype=object)

            _output[0][0] = "Sample ID"
            _output[0][1:] = self._list[0].masses

            for index, s in enumerate(self._list):
                _output[index + 1][0] = s.identifier
                _output[index + 1][1:] = s.intensities

            np.savetxt(fp, _output, delimiter=sep, fmt="%s")

        if output_type == "metaboanalyst":
            _to_metaboanalyst()
        elif output_type == "matrix":
            if self.value_imputated:
                _to_matrix()
            else:
                raise ValueError("SpectrumList must be Value Imputated first!")
        else:
            raise ValueError("%s not a valid output type" % (output_type))
