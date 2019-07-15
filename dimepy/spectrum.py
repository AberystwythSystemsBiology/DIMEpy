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
from typing import Tuple, List

import math
from scipy.stats import binned_statistic
from joblib import Parallel, delayed

import pandas as pd

from pymzml.run import Reader as pymzmlReader

from .utils import terms, bin_masses_and_intensities
from .scan import Scan


class Spectrum:
    def __init__(self,
                 filepath: str,
                 identifier: str = None,
                 injection_order: int = None,
                 stratification: str = None,
                 snr_estimator: str = False):
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

        self.read_scans = []
        self._masses = False
        self._intensities = False

        self._scans, self._to_use = self._base_load()

    def _base_load(self) -> Tuple[np.array, np.array]:
        extraAccessions = [[[y, ["value"]] for y in terms[x]]
                           for x in terms.keys()]

        # Flatten the list of lists of lists into a list of lists.
        extraAccessions = list(itertools.chain.from_iterable(extraAccessions))

        reader = pymzmlReader(self.filepath, extraAccessions=extraAccessions)

        scans = []
        to_use = []

        for scan in reader:
            scans.append(scan)
            to_use.append(True)

        return np.array(scans), np.array(to_use)

    def limit_polarity(self, polarity: str) -> None:
        def _determine_polarity(scan) -> str:
            scan_polarity = None
            for polarity_acc in terms["polarity"]:
                if scan.get(polarity_acc) != None:
                    scan_polarity = terms["polarity"][polarity_acc]

            return scan_polarity

        for index, scan in enumerate(self._scans):
            if _determine_polarity(scan) != polarity.upper():
                self._to_use[index] = False

    def limit_apex(self, mad_multiplyer: int = 1) -> None:

        tics = np.array([scan.TIC for scan in self._scans])
        mad = np.mean(np.absolute(np.array(tics) - np.mean(tics)))

        apex_index = tics >= mad * mad_multiplyer

        sel = np.where(self._to_use == True)[0]

        for indx, sel in enumerate(sel):
            self._to_use[sel] = apex_index[indx]

    def reset(self) -> None:
        self._to_use = self._to_use == False
        self._masses = False
        self._intensities = False
        self.read_scans = []

    def load_scans(self) -> None:
        scans = []
        masses = []
        intensities = []
        for scan in self._scans[self._to_use]:
            scans.append(Scan(scan, snr_estimator=self.snr_estimator))
        self.read_scans = scans

        self._get_masses_and_ints()

    def bin(self, bin_width: float = 0.01, statistic: str = "mean"):

        self._masses, self._intensities = bin_masses_and_intensities(
            self.masses, self.intensities, bin_width, statistic)

    def remove_spurious_peaks(self,
                              bin_width: float = 0.01,
                              threshold: float = 0.5,
                              scan_grouping: int = 50.0):
        def _determine_scan_group():
            medians = [np.mean(x.masses) for x in self.scans]

            bins = np.arange(np.min(medians), np.max(medians), scan_grouping)

            _, _, binnumber = binned_statistic(medians,
                                               medians,
                                               statistic="mean",
                                               bins=bins)

            scan_groups = {}

            for index, group in enumerate(np.unique(binnumber)):
                scan_groups[index] = self.scans[binnumber == group]

            return scan_groups

        def _get_bins(scan_list):
            mass_ranges = np.array([x.mass_range for x in scan_list])
            min_mass = np.min([x[0] for x in mass_ranges])
            max_mass = np.max([x[1] for x in mass_ranges]) + bin_width
            return np.arange(min_mass, max_mass, step=bin_width)

        def _calculate_bins(scan_list, bins):

            scan_index = {}

            for index, scan in enumerate(scan_list):

                _, _, binnumber = binned_statistic(scan.masses,
                                                   scan.intensities,
                                                   bins=bins)

                counts = []

                for bin_index in range(len(bins)):
                    number_bins = np.count_nonzero(binnumber == bin_index)
                    counts.append(number_bins)

                scan_index[index] = counts

            df = pd.DataFrame(scan_index).T
            df.columns = bins
            df = df.loc[:, (df != 0).any(axis=0)]

            df = df.replace(0, np.nan)

            to_keep = df.columns[df[df.columns].isnull().sum() <= len(df) *
                                 threshold]
            return to_keep

        def _remove_from_scans(scan_list, non_spurios_masses):
            for scan in scan_list:
                masses = []
                intensities = []
                for non_spurios_mass in non_spurios_masses:
                    keep = np.where(
                        np.logical_and(
                            scan.masses >= non_spurios_mass,
                            scan.masses <= non_spurios_mass + bin_width))
                    masses.extend(scan.masses[keep].tolist())
                    intensities.extend(scan.intensities[keep].tolist())
                scan.masses = masses
                scan.intensities = intensities

        self._masses = False
        self._intensities = False

        if scan_grouping:
            scan_groups = _determine_scan_group()
        else:
            scan_groups = {0: self.scans}

        for scan_group in scan_groups:
            scan_list = scan_groups[scan_group]
            bins = _get_bins(scan_list)
            non_spurios_masses = _calculate_bins(scan_list, bins)
            _remove_from_scans(scan_list, non_spurios_masses)

    def _get_masses_and_ints(self):
        masses = []
        intensities = []
        for scan in self.read_scans:
            masses.extend(scan.masses.tolist())
            intensities.extend(scan.intensities.tolist())

        masses = np.array(masses)
        intensities = np.array(intensities)

        sorted_idx = np.argsort(masses)

        self._masses = masses[sorted_idx]
        self._intensities = intensities[sorted_idx]

    @property
    def scans(self):
        return self.read_scans

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
            raise ValueError(
                "No intensities generated, run Spectrum.get first")

    @property
    def mass_range(self):
        return [np.min(self.masses), np.max(self.masses)]
