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
import pandas as pd
from .spectrum import Spectrum

class SpectrumList:

    def __init__(self):
        self._list = []

        self.binned = False


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

    def bin(self, bin_width: float = 0.5, statistic: str = "mean"):
        """
        Method to conduct mass binning to nominal mass and mass spectrum
        generation across a SpectrumList.

        Arguments:
            bin_width (float): The mass-to-ion bin-widths to use for binning.
            statistic (str): The statistic to use to calculate bin values.
        """
        def _get_mass_range():
            mass_range = [x.mass_range for x in self._list]
            min_mass = np.min(np.min([x[0] for x in mass_range])) - bin_width
            max_mass = np.max(
                np.max([x[1] for x in mass_range])) + bin_width
            return min_mass, max_mass

        def _calculate_bins(bins):

            bin_dict = {x: [] for x in bins}

            _intensities = []

            for index, spectrum in enumerate(self._list):
                masses = spectrum.masses
                intensities = spectrum.intensities

                binned_intensities, _, _ = binned_statistic(
                    masses, intensities, statistic=statistic, bins=bins)

                binned_masses, _, _ = binned_statistic(
                    masses, masses, statistic=statistic, bins=bins)

                index = ~np.isnan(binned_intensities)

                binned_masses = binned_masses[index]

                for bin_indx, bin in enumerate(bins[:-1][index]):
                    bin_dict[bin].append(binned_masses[bin_indx])

                _intensities.append(binned_intensities)
            
            return _intensities, bin_dict

        def calculate_masses(bin_dict):
            bins = []
            for bin in bin_dict:
                if bin_dict[bin] != []:
                    bins.append(np.mean(bin_dict[bin]))
                else:
                    bins.append(bin)
            return sorted(bins)

        min_mass, max_mass = _get_mass_range()

        bins = np.arange(min_mass, max_mass, step=bin_width)
        
        intensities, bin_dict = _calculate_bins(bins)
        masses = calculate_masses(bin_dict)
        
        for index, binned_ints in enumerate(intensities):
            print(binned_ints[0])

            not_null = np.where(binned_ints == np.nan)

            print(not_null)
