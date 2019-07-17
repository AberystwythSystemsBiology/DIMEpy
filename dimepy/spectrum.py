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

from scipy.stats import binned_statistic
from pymzml.run import Reader as pymzmlReader

from .scan import Scan
from .utils import terms, bin_masses_and_intensities


class Spectrum:

    def __init__(self,
                 filepath: str,
                 identifier: str,
                 injection_order: int = None,
                 stratification: str = None,
                 snr_estimator: str = False):
        """
        Initialise Spectrum object for a given mzML file.

        Arguments:
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
        extraAccessions = [
            [[y, ["value"]] for y in terms[x]] for x in terms.keys()
        ]

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
        """
        Limit the Scans found within the mzML file to whatever polarity is given.
        This is quite useful for those of us that are using mixed polarity experiments.

        Arguments:
            polarity (string): polarity type of the scans required (positive/negative)
        """

        def _determine_polarity(scan) -> str:
            scan_polarity = None
            for polarity_acc in terms["polarity"]:
                if scan.get(polarity_acc) != None:
                    scan_polarity = terms["polarity"][polarity_acc]

            return scan_polarity

        for index, scan in enumerate(self._scans):
            if _determine_polarity(scan) != polarity.upper():
                self._to_use[index] = False

    def limit_infusion(self, mad_multiplyer: int = 1) -> None:
        """
        This method is a slight extension of the work by Manfred Beckmann
        (meb@aber.ac.uk) in FIEMSpro in which we use the mean absolute
        deviation to determine when the infusion has taken place.

        Infusion Profile (Sketch):

               _
              / \
             /   \_
        ____/       \_________________
        0     0.5     1     1.5     2 [min]
            |--------| Apex

        
        Arguments:
            mad_multiplier (int): The multiplier for the mean absolute
            deviation method to take the infusion profile from.
        
        """
        tics = np.array([scan.TIC for scan in self._scans])
        mad = np.mean(np.absolute(np.array(tics) - np.mean(tics)))

        apex_index = tics >= mad * mad_multiplyer

        sel = np.where(self._to_use == True)[0]

        for indx, sel in enumerate(sel):
            self._to_use[sel] = apex_index[indx]

    def reset(self) -> None:
        """
        A method to reset the Spectrum object completely.        
        """
        self._to_use = self._to_use == False
        self._masses = False
        self._intensities = False
        self.read_scans = []

    def load_scans(self) -> None:
        """
        A method to load the scans in.

        Note: 
            If you are using mixed polarity methods, or only require the
            infusion profile - please ensure you run limit_infusion and
            limit_profile.
        """
        scans = []

        for scan in self._scans[self._to_use]:
            scan = Scan(scan, snr_estimator=self.snr_estimator)
            scans.append(scan)

        self.read_scans = scans
        self._load_masses_and_ints_from_scans()

    def _load_masses_and_ints_from_scans(self) -> None:

        masses = []
        intensities = []

        for scan in self.read_scans:
            masses.extend(scan.masses)
            intensities.extend(scan.intensities)

        masses = np.array(masses)
        intensities = np.array(intensities)

        sorted_idx = np.argsort(masses)

        self._masses = masses[sorted_idx]
        self._intensities = intensities[sorted_idx]

    def bin(self, bin_width: float = 0.01, statistic: str = "mean"):
        """
        Method to conduct mass binning to nominal mass and mass spectrum
        generation.

        Arguments:
            bin_width (float): The mass-to-ion bin-widths to use for binning.
            statistic (str): The statistic to use to calculate bin values.
        """
        self._masses, self._intensities = bin_masses_and_intensities(
            self.masses, self.intensities, bin_width, statistic)

    def remove_spurious_peaks(self,
                              bin_width: float = 0.01,
                              threshold: float = 0.25,
                              scan_grouping: float = 50.0):
        """
        Method that's highly influenced by Jasen Finch's (jsf9@aber.ac.uk)
        binneR, in which spurios peaks can be removed. At the time of writing,
        this method has serious performance issues and needs to be rectified.

        Arguments:
            bin_width (float): The mass-to-ion bin-widths to use for binning.
            threshold (float): Percentage of scans in which a peak must be in
            in order for it to be considered.
            scan_grouping (float): Mass-to-ion scan groups.


        Note:
            load_scans() must first be run in order for this to work.

        """

        def _determine_scan_group():
            medians = [np.mean(x.masses) for x in self.scans]

            bins = np.arange(np.min(medians), np.max(medians), scan_grouping)

            _, _, binnumber = binned_statistic(medians,
                                               medians,
                                               statistic="mean",
                                               bins=bins)

            scan_groups = {}

            for index, group in enumerate(np.unique(binnumber)):
                scan_groups[index] = np.array(self.scans)[binnumber == group]

            return scan_groups

        def _get_bins(scan_list: list) -> np.array:
            mass_ranges = np.array([x.mass_range for x in scan_list])
            min_mass = np.min([x[0] for x in mass_ranges]) - bin_width
            max_mass = np.max([x[1] for x in mass_ranges]) + bin_width
            return np.arange(min_mass, max_mass, step=bin_width)

        def _calculate_bins(scan_list: list, bins: np.array) -> np.array:

            scan_index = []

            for index, scan in enumerate(scan_list):

                _, _, binnumber = binned_statistic(scan.masses,
                                                   scan.intensities,
                                                   bins=bins)

                counts = []

                for bin_index in range(len(bins)):
                    number_bins = np.count_nonzero(binnumber == bin_index)
                    counts.append(number_bins)

                scan_index.append(np.array(counts))

            _tmp_si = np.array(scan_index)
            return bins[_tmp_si.sum(axis=0) >= len(scan_list) * threshold]

        def _remove_from_scans(scan_list: list,
                               non_spurios_masses: np.array) -> None:
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

        if scan_grouping:
            scan_groups = _determine_scan_group()
        else:
            scan_groups = {0: self.scans}

        for scan_group in scan_groups:
            scan_list = scan_groups[scan_group]
            bins = _get_bins(scan_list)
            non_spurios_masses = _calculate_bins(scan_list, bins)
            _remove_from_scans(scan_list, non_spurios_masses)

        # Load in new masses and intensities.
        self._load_masses_and_ints_from_scans()

    @property
    def scans(self) -> List[Scan]:
        return self.read_scans

    @property
    def masses(self) -> np.array:
        if type(self._masses) != bool:
            return self._masses
        else:
            raise ValueError("No masses generated, run Spectrum.get first.")

    @property
    def intensities(self) -> np.array:
        if type(self._intensities) != bool:
            return self._intensities
        else:
            raise ValueError("No intensities generated, run Spectrum.get first")

    @property
    def mass_range(self) -> Tuple[float, float]:
        return [np.min(self.masses), np.max(self.masses)]
