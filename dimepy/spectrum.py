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
            if index > 10:
                break
            scan = Scan(pymzmlSpectrumInstance, snr_estimator=self.snr_estimator)

            scans.append(scan)
            to_use.append(True)

        return np.array(scans), np.array(to_use)

            

    def get_polarity(self, polarity: str):
        for index, scan in enumerate(self._scans):
            if scan.polarity != polarity.upper():
                self._to_use[index] = False

    @property
    def scans(self):
        return self._scans[self._to_use == True]