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
import math
from scipy.stats import binned_statistic

import pandas as pd

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

            if index > 5:
                break

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

    def bin(self, bin_width: float = 0.01, statistic: str = "mean"):
        min_mass, max_mass = self.mass_range

        min_mass = math.floor(min_mass)
        max_mass = math.ceil(max_mass)+bin_width

        bins = np.arange(min_mass, max_mass, bin_width)

        binned_intensities, _, _ = binned_statistic(
            self.masses,
            self.intensities,
            statistic = statistic,
            bins = bins
        )

        binned_masses, _, _= binned_statistic(
            self.masses,
            self.masses,
            statistic = statistic,
            bins = bins
        )

        index = ~np.isnan(binned_intensities)

        self._masses = binned_masses[index]
        self._intensities = binned_intensities[index]

    def _calculate_mass_range(self):
        return [np.min(self.masses), np.max(self.masses)]


    def remove_spurious_peaks(self, bin_width: float = 0.01, threshold: float = 0.5, scan_grouping: int = 50.0):
        
        def _determine_scan_group():
            medians = [np.mean(x.masses) for x in self.scans]

            bins = np.arange(np.min(medians), np.max(medians), scan_grouping)

            _, _, binnumber = binned_statistic(
                medians,
                medians,
                statistic = "mean",
                bins = bins
            )

            scan_groups = {}

            for index, group in enumerate(np.unique(binnumber)):
                scan_groups[index] = self.scans[binnumber == group]

            return scan_groups


        def _get_bins(scan_list):
            mass_ranges = np.array([x.mass_range for x in scan_list])
            min_mass = np.min([x[0] for x in mass_ranges])
            max_mass = np.max([x[1] for x in mass_ranges])+bin_width
            return np.arange(min_mass, max_mass, step=bin_width)


        def _calculate_bins(scan_list, bins):

            scan_index = {}

            for index, scan in enumerate(scan_list):

                _, _, binnumber = binned_statistic(
                    scan.masses,
                    scan.intensities,
                    bins = bins
                )

                counts = []

                for bin_index in range(len(bins)):
                    number_bins = np.count_nonzero(binnumber == bin_index)
                    counts.append(number_bins)


                scan_index[index] = counts

            df = pd.DataFrame(scan_index).T
            df.columns = bins
            df = df.loc[:, (df != 0).any(axis=0)]

            df = df.replace(0, np.nan)

            to_keep = df.columns[df[df.columns].isnull().sum() <= len(df)*threshold]
            return to_keep
        
        def _remove_from_scans(scan_list, non_spurios_masses):
            for scan in scan_list:
                masses = []
                intensities = []
                for non_spurios_mass in non_spurios_masses:
                    keep = np.where(np.logical_and(scan.masses>=non_spurios_mass, scan.masses<=non_spurios_mass+bin_width))
                    masses.extend(scan.masses[keep].tolist())
                    intensities.extend(scan.intensities[keep].tolist())
                scan.masses = masses
                scan.intensities = intensities

        self._masses = False
        self._intensities = False

        if scan_grouping:
            scan_groups = _determine_scan_group()
        else:
            scan_groups = {0 : self.scans}

        for scan_group in scan_groups:
            scan_list = scan_groups[scan_group]
            bins = _get_bins(scan_list)
            non_spurios_masses = _calculate_bins(scan_list, bins)
            _remove_from_scans(scan_list, non_spurios_masses)




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
            raise ValueError("No intensities generated, run Spectrum.get first")

    @property
    def mass_range(self):
        return self._calculate_mass_range()
