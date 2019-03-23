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

import itertools
import numpy as np
import matplotlib.pyplot as plt

from pymzml.run import Reader as pymzmlReader

from .utils import terms
from .scan import Scan


class Spectrum:

    def __init__(self, filepath: str, identifier: str = None, injection_order: int = None, stratification: str = None, snr_estimator: str = False):
        """
        Initialise Spectrum object for a given mzML file.

        Args:
            filepath (str): Path to the mzML file to parse.
            identifier (str): Unique identifier for the Spectrum object.
            injection_order (int): The injection number of the Spectrum object.
            stratification (str): Class label of the Spectrum object)
            snr_estimator (str): Signal to noise method used to filter.
        """
        self.filepath = filepath
        self.identifier = identifier
        self.injection_order = injection_order
        self.stratification = stratification
        self.snr_estimator = snr_estimator

        self._scans, self._to_use = self.load()

        self._masses = False
        self._intensities = False


    def load(self):
        extraAccessions=[
            [[y, ["value"]] for y in terms[x]] for x in terms.keys()
        ]

        # Flatten the list of lists of lists into a list of lists.
        extraAccessions = list(itertools.chain.from_iterable(extraAccessions))

        reader = pymzmlReader(self.filepath, extraAccessions=extraAccessions)

        scans = []
        to_use = []

        for index, pymzmlSpectrumInstance in enumerate(reader):
            scan = Scan(pymzmlSpectrumInstance, snr_estimator=self.snr_estimator)

            scans.append(scan)
            to_use.append(True)

        return np.array(scans), np.array(to_use)


    def get_polarity(self, polarity: str):
        for index, scan in enumerate(self._scans):
            if scan.polarity != polarity.upper():
                self._to_use[index] = False


    def get_apex(self, mad_multiplyer: int = 1, filename: str = False):
        tics = np.array([scan.total_ion_count for scan in self.scans])
        mad = np.mean(np.absolute(np.array(tics) - np.mean(tics)))

        apex_index = tics >= mad * mad_multiplyer

        if filename:
            plt.figure()
            plt.title("Apex Plot: %s" % self.identifier)
            plt.plot(tics)
            plt.ylabel("Total Ion Count")
            plt.xlabel("Scan Number")
            plt.tight_layout()
            plt.savefig(filename)

        sel = np.where(self._to_use == True)[0]

        for indx, sel in enumerate(sel):
            self._to_use[sel] = apex_index[indx]


    def reset(self):
        self._to_use = self._to_use[np.where(self._to_use == False)] == True
        self._masses = False
        self._intensities = False

    def get(self):
        masses = []
        intensities = []

        for scan in self.scans:
            masses.extend(scan.masses)
            intensities.extend(scan.intensities)

        spectrum = list(zip(masses, intensities))

        # Sort by masses
        spectrum.sort(key=lambda x: float(x[0]))

        self._masses = np.array([x[0] for x in spectrum])
        self._intensities = np.array([x[1] for x in spectrum])

    @property
    def scans(self):
        return self._scans[self._to_use == True]

    @property
    def masses(self):
        if type(self._masses) != bool:
            return self._masses
        else:
            raise ValueError("No masses generated, run Spectrum.get first.")

    @property
    def intensities(self):
        if type(self._intensities) != bool:
            return self._intensities
        else:
            raise ValueError("No intensities generated, run get Spectrum.first.")
