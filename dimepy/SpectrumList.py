# -*- coding: utf-8 -*-
import pandas as pd
import cPickle as pkl
import collections
import csv
from copy import copy
import numpy as np
from Spectrum import Spectrum
import matplotlib.pyplot as plt

class SpectrumList(object):
    """An object to hold and transform spectrum objects.

    Parameters
    ----------

    _spectrum : list, optional (default=[])

    """
    _outlier_detected = False
    _binned = False
    _value_imputated = False
    _scaled = False
    _value_imputated = False

    def __init__(self, _spectrum=[]):
        self._spectrum = _spectrum

    def outlier_detection(self, threshold=3, inplace=True, plot=False):
        """Perform outlier detection.

        This method detects and removes sample outliers through the use of the
        mean absolute deviation over the total ion counts.

        Parameters
        ----------

        threshold : float, optional (default=3)

        plot : boolean, optional (default=False)

        inplace : boolean, optional (default=True)
            If False then return the outlier detected SpectrmList, else make the
            change within the SpectrmList.

        """

        def _calculate_mad():
            tics = np.array([s.tic for s in self._spectrum])
            return np.mean(np.absolute(tics - np.mean(tics)))

        def _find_outliers(mad):
            return [s for s in self._spectrum if s.tic >= mad*threshold]

        def _plot(outliers):
            fig, ax = plt.subplots()
            for outlier in outliers:
                plt.scatter(outlier.injection_order, outlier.tic, color="red")
                ax.annotate(outlier.id, (outlier.injection_order, outlier.tic))
            for spectrum in self._spectrum:
                if spectrum not in outliers:
                    plt.scatter(spectrum.injection_order, spectrum.tic, color="blue")
                    ax.annotate(spectrum.id, (spectrum.injection_order, spectrum.tic))
            for t in range(threshold):
                l = "%s times MAD" % (t+1)
                ax.axhline(mad*(t+1), linestyle="--", color="rgboy"[t], label=l)
            plt.ylabel("Injection Order")
            plt.xlabel("Total Ion Count")
            plt.tight_layout()
            plt.legend(loc="best")
            plt.show()

        mad = _calculate_mad()
        outliers = _find_outliers(mad)

        if plot == True:
            _plot(outliers)

        if inplace == True:
            self.delete(outliers)
        else:
            return self.delete(outliers, inplace=False)

    def append(self, s):
        """Append Spectrum to the end of the SpectrumList.

        Parameters
        ---------

        s : Spectrum
            A Spectrum object to append to the SpectrumList.

        """
        if type(s) == Spectrum:
            self._spectrum.append(s)
        else:
            raise ValueError("%s is not a valid object type, need Spectrum object" % type(s))

    def to_pickle(self, fp):
        """Dump the SpectrmList to a pickled object.

        Parameters
        ----------

        fp : string
            File path to write the SpectrmList to.
        """
        with open(fp, "wb") as outfile:
            pkl.dump(self, outfile, pkl.HIGHEST_PROTOCOL)

    def delete(self, s, inplace=True):
        """Remove a given Spectrum object from the list.

        Parameters
        ----------

        s : Spectrum
            A spectrum object or a list of spectrum objects to be removed.

        inplace : boolean, optional (default=True)
            If False then return the corrected intensities array, else make the
            change within the object.

        """
        if inplace == True:
            if type(s) == Spectrum:
                self._spectrum.remove(s)
            elif type(s) == type([]):
                [self._spectrum.remove(x) for x in s]

        else:
            cs = copy(self._spectrum)
            if type(s) == Spectrum:
                cs.remove(s)
            elif type(s) == type([]):
                [cs.remove(x) for x in s]
            return cs

    def to_csv(self, fp, delimiter=","):
        """Write the SpectrmList to a delimited file.

        Parameters
        ----------

        fp : string
            File path to write the delimited file to.

        delimiter : string, optional (default=",")
            Field delimiter for the output file.

        """
        raw_values = []

        for spectrum in self._spectrum:
            raw_values.append([spectrum.id] + spectrum.masses.tolist())
            raw_values.append([spectrum.label] + spectrum.intensities.tolist())
        raw_values = pd.DataFrame(raw_values).T.replace(
            np.nan, "", regex=True).values

        with open(fp, "wb") as output:
            writer = csv.writer(output, delimiter=delimiter)
            for row in raw_values:
                writer.writerow(row)
            output.close()

    def flatten_to_dataframe(self):
        output = []
        for spectrum in self.tolist():
            df = pd.DataFrame(spectrum.intensities).T
            df.columns = spectrum.masses
            df.index = [spectrum.id]
            output.append(df)
        return pd.concat(output, axis=0)

    def tolist(self):
        return self._spectrum

    def pandas_to_spectrum(self, df):
        masses = df.columns
        for id, values in df.iterrows():
            intensities = values.values
            spectrum = [x for x in self.to_list() if x.id == id][0]
            spectrum.masses = masses
            spectrum.intensities = intensities

    def get_mass_range(self):
        smallest = None
        largest = None

        for spectrum in self.tolist():
            if min(spectrum.masses) < smallest or smallest is None:
                smallest = min(spectrum.masses)
            if max(spectrum.masses) > largest or largest is None:
                largest = max(spectrum.masses)

        return smallest, largest

    def __repr__(self):
        return repr(self.to_list())
