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


from .utils import terms
import numpy as np

class Scan:
    def __init__(self,
                 pymzml_spectrum,
                 snr_estimator: str = False,
                 peak_type: str = "raw"):
        """
            Initalise a Scan object for a given pymzML Spectrum.

            Args:
                pymzml_spectrum (pymzml.Spectrum): Spectrum object

        """
        self.pymzml_spectrum = pymzml_spectrum
        self.polarity = self._get_polarity()

        self.peak_type = peak_type

        if snr_estimator:
            self.snr = self._estimate_snr(snr_estimator)
        else:
            self.snr = False

        self.masses, self.intensities = self._get_spectrum()

        self.total_ion_count = np.sum(self.intensities)
        
        self.mass_range = np.array(
            np.min(self.masses),
            np.max(self.masses)
        )

    def _estimate_snr(self, snr_estimator):
        return self.pymzml_spectrum.estimated_noise_level(mode=snr_estimator)


    def _get_spectrum(self):
        try:
            peaks = getattr(self.pymzml_spectrum, "peaks")(self.peak_type)
            
            masses, intensities = [np.array(x) for x in zip(*peaks)]

            if self.snr:
                not_noise = intensities >= self.snr
                masses = masses[not_noise]
                intensities = intensities[not_noise]

            return np.array(masses), np.array(intensities)

        except ValueError:
            raise ValueError("%s is not a supported peak type." % (self.peak_type))

    def _get_polarity(self):
        polarity = None
        for polarity_accession in terms["polarity"].keys():
            if self.pymzml_spectrum.get(polarity_accession) != None:
                polarity = terms["polarity"][polarity_accession]
        return polarity