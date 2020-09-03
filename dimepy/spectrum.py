# -*- coding: utf-8 -*-
# encoding: utf-8

# Copyright (c) 2017-2020 Keiron O'Shea
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
import logging
from typing import Tuple, List
from scipy.stats import binned_statistic
from pymzml.run import Reader as pymzmlReader
from .scan import Scan
from .utils import terms, bin_masses_and_intensities
import itertools
import os
import matplotlib.pyplot as plt

"""
The class :py:class:`Spectrum` has been designed to load data
from a mzML file.

.. note::
    This class is still being actively developed and will likely change over time.
"""

class Spectrum(object):

    """ Initialise Spectrum object for a given mzML file.

    """

    def __init__(self,
                 filepath: str,
                 identifier: str = None,
                 injection_order: int = None,
                 stratification: str = None,
                 snr_estimator: str = False,
                 peak_type: str = "raw",
                 is_qc: bool= False,
                 MS1_precision: float = 5e-6,
                 MSn_precision: float = 20e-6
                 ):

        """

        Initialise a Spectrum object for a given mzML file.

        Arguments:
            filepath (str): Path to the mzML file to parse.

            identifier (str): Unique identifier for the Spectrum object. If none is given,
                it will default to the filepath.

            injection_order (int): The injection number of the Spectrum object.

            stratification (str): Class label of the Spectrum object.

            snr_estimator (str): Signal to noise method used to filter.

                Currently supported signal-to-noise estimation methods are:
                    * 'median' (default)
                    * 'mean'
                    * 'mad'

            peak_type (str): What peak type to load in.

                Currently supported peak types are:
                    * raw (default)
                    * centroided
                    * reprofiled

            is_qc (bool): Whether the Sample is QC.

            MS1_precision (float): Measured precision for the MS level 1.

            MSn_precision (float): Measured precision for the MS level n.

    """
        self.filepath = filepath
        self.identifier = identifier

        if identifier == None:
            self.identifier = os.path.basename(filepath)

        self.injection_order = injection_order
        self.stratification = stratification
        self.snr_estimator = snr_estimator
        self.peak_type = peak_type
        self.is_qc = is_qc
        self.MS1_precision = MS1_precision
        self.MSn_precision = MSn_precision

        self.read_scans = []
        self._masses = False
        self._intensities = False

        self._scans, self._to_use = self._base_load()

    def _base_load(self) -> Tuple[np.array, np.array]:

        logging.info("Loading %s" % (self.filepath))

        extraAccessions = [
            [[y, ["value"]] for y in terms[x]] for x in terms.keys()
        ]

        # Flatten the list of lists of lists into a list of lists.
        extraAccessions = list(itertools.chain.from_iterable(extraAccessions))

        reader = pymzmlReader(
            self.filepath,
            extraAccessions=extraAccessions,
            MS1_Precision = self.MS1_precision,
            MSn_Precision = self.MSn_precision
        )

        scans = []
        to_use = []

        for scan in reader:
            scans.append(scan)
            to_use.append(True)

        return np.array(scans), np.array(to_use)

    def limit_polarity(self, polarity: str, verbose: bool = False) -> None:
        """
        Limit the Scans found within the mzML file to whatever polarity is given.
        This should only be called where fast-polarity switching is used.

        Arguments:
            polarity (str): polarity type of the scans required

                Supported polarity types are:
                    * 'positive'
                    * 'negative'

            verbose (bool): enable verbose output.

        """

        if polarity.upper() not in ["POSITIVE", "NEGATIVE"]:
            raise AttributeError("%s not a valid option" % (polarity))

        def _determine_polarity(scan) -> str:
            scan_polarity = None
            for polarity_acc in terms["polarity"]:
                if scan.get(polarity_acc) != None:
                    scan_polarity = terms["polarity"][polarity_acc]
            return scan_polarity

        for index, scan in enumerate(self._scans):
            if _determine_polarity(scan) != polarity.upper():
                self._to_use[index] = False
                logging.info("Scan %i is not %s polarity" % (index, polarity))

            if verbose and self._to_use[index]:
                print("Scan %i is %s polarity" % (index, polarity))


    def limit_infusion(self, threshold: int = 3, plot = False) -> None:
        """
        This method is a slight extension of Manfred Beckmann's (meb@aber.ac.uk)
        LCT/Q-ToF scan retrieval method in FIEMSpro in which we use the median absolute
        deviation of all TICs within a Spectrum to determine when
        the infusion has taken place.

        Consider the following Infusion Profile:
        ::
                 _
                / \ 
               /   \_
          ____/       \_________________
          0     0.5     1     1.5     n [scan number]
              |--------| Apex
        
        We are only interested in the scans in which the infusion takes place
        (20 - 50 seconds). Applying this method changes the to_use values to only
        be True where the TIC is >= TIC * mad_multiplier.
        
        Arguments:

            mad_multiplier (int): The multiplier for the median absolute
                deviation method to take the infusion profile from.

            plot (bool): Plot the results.

        """

        def _calculate_mad(tics: np.array) -> float:
            return np.median(np.abs(tics - np.median(tics)))

        def _get_mask(tics: np.array, mad: float) -> np.array:
            tics = tics[:, None]
            median = np.median(tics, axis=0)
            diff = np.sum((tics - median)**2, axis=-1)
            diff = np.sqrt(diff)

            med_abs_deviation = np.median(diff)

            modified_z_score = 0.6745 * diff / med_abs_deviation

            return modified_z_score >= threshold

        def _plot(tics: np.array, mad: float, ini_scan: int, end_scan: int):
            plt.figure()
            plt.title("Apex Selection Plot \n%s" % (self.identifier))
            plt.plot(tics)
            plt.xlim(0, len(tics))
            plt.ylim(0, max(tics * 1.1))
            plt.axvline(ini_scan, ls="--", c="red")
            plt.axvline(end_scan, ls="--", c="red")
            plt.xlabel("Scan Number")
            plt.ylabel("Total Ion Current (TIC)")
            plt.tight_layout()
            if type(plot) == bool:
                plt.show()
            else:
                plt.savefig(plot)

        tics = np.array([scan.TIC for scan in self._scans[self._to_use]])
        mad = _calculate_mad(tics)
        apex_index = _get_mask(tics, mad)

        logging.info("Mean Absolute Deviation Calculated as %f TIC" % (mad))

        ini_scan = np.min(np.where(apex_index == True)[0])
        end_scan = np.max(np.where(apex_index == True)[0])

        if plot:
            _plot(tics, mad, ini_scan, end_scan)

        to_use = np.where(self._to_use == True)[0]

        for i, j in enumerate(to_use):
            self._to_use[j] = apex_index[i]

    def reset(self) -> None:
        """
        A method to reset the Spectrum object in its entirety.
        """
        self._to_use = self._to_use == False
        self._masses = False
        self._intensities = False
        self.read_scans = []

    def load_scans(self) -> None:
        """
        This method loads the scans in accordance to whatever Scans are
        set to True in the to_use list.

        .. note:: If you want to actually make use of masses and intensities
            (you probably do), then ensure that you call this method.
        """
        scans = []

        for scan in self._scans[self._to_use]:
            scan = Scan(scan, snr_estimator=self.snr_estimator, peak_type=self.peak_type)
            scans.append(scan)

        self.read_scans = scans
        self._load_masses_and_ints_from_scans()

    def _load_masses_and_ints_from_scans(self) -> None:

        masses = []
        intensities = []

        for scan in self.read_scans:
            for m, i in zip(scan.masses, scan.intensities):
                if i > 0.0:
                    masses.append(m)
                    intensities.append(i)

        masses = np.array(masses)
        intensities = np.array(intensities)

        sorted_idx = np.argsort(masses)

        self._masses = masses[sorted_idx]
        self._intensities = intensities[sorted_idx]

    def bin(self, bin_width: float = 0.01, statistic: str = "mean"):
        """"
        Method to conduct mass binning to nominal mass and mass spectrum
        generation across a Spectrum.

        Arguments:
            bin_width (float): The mass-to-ion bin-widths to use for binning.

            statistic (str): The statistic to use to calculate bin values.

                Supported statistic types are:
                    * 'mean' (default): compute the mean of intensities for points within each bin.
                        Empty bins will be represented by NaN.
                    * 'std': compute the standard deviation within each bin. This is
                        implicitly calculated with ddof=0.
                    * 'median': compute the median of values for points within each bin.
                        Empty bins will be represented by NaN.
                    * 'count': compute the count of points within each bin.
                        This is identical to an unweighted histogram. values array is not referenced.
                    * 'sum': compute the sum of values for points within each bin.
                        This is identical to a weighted histogram.
                    * 'min': compute the minimum of values for points within each bin.
                        Empty bins will be represented by NaN.
                    * 'max': compute the maximum of values for point within each bin.
                        Empty bins will be represented by NaN.
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
        but should still work as intended (provided that you don't mind how long
        it takes to complete)

        Arguments:
            bin_width (float): The mass-to-ion bin-widths to use for binning.

            threshold (float): Percentage of scans in which a peak must be in
                in order for it to be considered.

            scan_grouping (float): Mass-to-ion scan groups, this splits the
                scans into groups to ease the processing somewhat. It
                is strongly recommended that you keep this at it's default
                value of of 50.0

        .. note:: load_scans() must first be run in order for this to work.
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
            raise ValueError("No masses generated, run load_scans() first!")

    @property
    def intensities(self) -> np.array:
        if type(self._intensities) != bool:
            return self._intensities
        else:
            raise ValueError("No intensities generated, run load_scans() first!")
    
    @property
    def TIC(self) -> float:
        return sum(self.intensities)

    @property
    def to_use(self) -> List[bool]:
        return self._to_use

    @property
    def mass_range(self) -> Tuple[float, float]:
        return [np.min(self.masses), np.max(self.masses)]
