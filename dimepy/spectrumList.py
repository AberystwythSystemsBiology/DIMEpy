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


class SpectrumList:

    def __init__(self):
        self._list = []

        self.binned = False
        self.normalised = False
        self.transformed = False


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

        pass

    def normalise(self, method: str = "tic") -> None:
        """
        Method to conduct sample independent intensity normalisation.

        Arguments:
            method (str): The normalisation method to use.
        """
        def _normie(spec: Spectrum):
            i = spec.intensities

            if method.upper() == "TIC":
                spec._intensities = np.divide(i,  np.sum(i)) * 1000
            elif method.upper() == "MEDIAN":
                spec._intensities = i - np.median(i) * 1000
            elif method.upper() == "MEAN":
                spec._intensities = i - np.mean(i) * 1000
            else:
                raise ValueError("%s is not a valid normalisation method" % (method))

        if self.normalised:
            for spec in self._list:
                _normie(spec)

            self.normalised = True
        else:
            raise ValueError("It looks like you've already normalised this data.")

    def transform(str, method: str = "log10") -> None:
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
                spec._intensities = np.log2(i + np.sqrt(i ** 2 + m ** 2)) / 2
            elif method.upper() == "SQRTasinh":
                spec._intensities = np.sqrt(i)
            elif method.upper() == "IHS":
                spec._intensities = np.array([math.asinh(x) for x in i])

        if self.transformed:
            for spec in self._list:
                _transform(spec)
            self.transformed = True
        else:
            raise ValueError("It looks like you've already transformed this data.")