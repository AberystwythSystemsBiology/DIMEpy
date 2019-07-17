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

from typing import Tuple
from scipy.stats import binned_statistic
from scipy.ndimage import find_objects
import numpy as np

terms = {"polarity": {"MS:1000129": "NEGATIVE", "MS:1000130": "POSITIVE"}}


def bin_masses_and_intensities(masses: np.array,
                               intensities: np.array,
                               bin_width: float = 0.25,
                               statistic: str = "mean",
                               bins: list = None) -> Tuple[np.array, np.array]:
    if bins == None:
        bins = np.arange(
            np.min(masses) - bin_width,
            np.max(masses) + bin_width, bin_width)

    statistic, _, bin_number = binned_statistic(masses,
                                                intensities,
                                                statistic=statistic,
                                                bins=bins)

    bin_number = bin_number

    masses = [np.mean(masses[x]) for x in find_objects(bin_number) if x != None]
    return np.array(masses), statistic[~np.isnan(statistic)]
