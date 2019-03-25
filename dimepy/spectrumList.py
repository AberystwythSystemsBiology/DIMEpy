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

class SpectrumList:
    def __init__(self, spectrum_list: list = [], bin_width=0.01, statistic="mean"):
        self.bin_width = bin_width
        self.statistic = statistic
        if len(spectrum_list) >= 2:
            self._spectrum_list = spectrum_list
            bin()
        else:
            # TODO: Throw exception if empty list passed
            exit(0)


    def bin(self):
        def _get_mass_range():
            mass_range = [x.mass_range for x in self._spectrum_list]
            min_mass = math.floor(np.min([x[0] for x in mass_range]))
            max_mass = math.ceil(np.max([x[1] for x in mass_range]) + self.bin_width)
            return min_mass, max_mass

        def _calculate_binning(bins):
            for spectrum in self._spectrum_list:
                masses = spectrum.masses
                intensities = spectrum.intensities

                binned_intensities, _, _ = binned_statistic(
                    masses,
                    intensities,
                    statistic=self.statistic,
                    bins=bins
                    )


        min_mass, max_mass = _get_mass_range()
        bins = np.arange(min_mass, max_mass, step=self.bin_width)


    def tolist(self):
        return self._spectrum_list
