# -*- coding: utf-8 -*-

from os.path import splitext
import pymzml
import numpy as np


class Scans(object):
    """A Scans object, containing scan information for a given file.

    Note
    ----

    Not supposed to be implemented outside of DIMEpy.

    """

    def __init__(self, fp, snr_estimator, max_snr, type):
        self.fp = fp
        self.snr_estimator = snr_estimator
        self.max_snr = max_snr
        self.type = type
        self.__scans = []
        self.__tics = []
        self.__polarities = []
        self._run()

    def _run(self):
        """

        """
        if splitext(self.fp)[1].upper() == ".MZML":
            self.__from_mzml()
        else:
            raise NotImplementedError(
                "DIMEpy currently only supports mzML files")

    def __from_mzml(self):
        """

        """

        def __gen_reader(eA):
            return pymzml.run.Reader(self.fp, extraAccessions=eA)

        def __get_ints_and_masses(scan):
            try:
                masses, intensities = zip(*getattr(scan, self.type))
            except AttributeError:
                raise ValueError("%s is not a supported peak type" % self.type)
            return np.array(masses), np.array(intensities)

        def __which_polarity(scan):
            polarity = None
            for pol, _ in eA:
                if scan.get(pol) is not None:
                    if pol == "MS:1000129":
                        polarity = "NEGATIVE"
                    elif pol == "MS:1000130":
                        polarity = "POSITIVE"
            return polarity

        def __apply_snr_filtering(scan, masses, intensities):
            indx = np.divide(
                intensities, scan.estimatedNoiseLevel(mode=self.snr_estimator))
            gq = indx >= self.max_snr
            masses = masses[gq]
            intensities = intensities[gq]
            return masses, intensities

        def __get_scans(reader, eA):
            for scan in reader:
                if scan["ms level"] is not None:
                    masses, intensities = __get_ints_and_masses(scan)
                    if self.snr_estimator is not None:
                        masses, intensities = __apply_snr_filtering(
                            scan, masses, intensities)
                    polarity = __which_polarity(scan)
                    self.__scans.append([masses, intensities])
                    self.__tics.append(sum(intensities))
                    self.__polarities.append(polarity)

        eA = [["MS:1000129", ["value"]], ["MS:1000130", ["value"]]]
        reader = __gen_reader(eA)
        __get_scans(reader, eA)

    def __from_raw(self):
        """

        """
        raise NotImplementedError("This has yet to be implemented.")

    def limiter(self, indx):
        """

        """
        self.__scans = np.array(self.__scans)[indx].tolist()
        self.__polarities = np.array(self.__polarities)[indx].tolist()
        self.__tics = np.array(self.__tics)[indx].tolist()

    @property
    def scans(self):
        return np.array(self.__scans)

    @property
    def polarities(self):
        return np.array(self.__polarities)

    @property
    def tics(self):
        return np.array(self.__tics)
