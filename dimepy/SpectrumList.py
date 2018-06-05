# -*- coding: utf-8 -*-
import pandas as pd
import cPickle as pkl
import collections
import csv
from scipy.stats import binned_statistic
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

    def __init__(self, _spectrum=[], dir=None):
        self._spectrum = _spectrum
        self._outlier_detected = False
        self._binned = False
        self._normalised = False
        self._transformed = False
        self._value_imputated = False
        self._scaled = False
        self._value_imputated = False
        if dir is not None:
            pass

    def outlier_detection(self, threshold=3, inplace=True, plot=False):
        """Perform outlier detection.

        This method detects and removes sample outliers through the use of the
        mean absolute deviation over the total ion counts.

        Parameters
        ----------

        threshold : float, optional (default=3)
            Value to calculate outlier threshold (the mean absolute deviation
            times the given threshold)

        plot : boolean, optional (default=False)
            Whether or not to plot the outlier detection process.

        inplace : boolean, optional (default=True)
            If False then return the outlier detected SpectrmList, else make
            the change within the SpectrmList.

        """

        def _calculate_mad():
            tics = np.array([s.tic for s in self._spectrum])
            return np.mean(np.absolute(tics - np.mean(tics)))

        def _find_outliers(mad):
            return [s for s in self._spectrum if s.tic >= mad * threshold]

        def _plot(outliers):
            fig, ax = plt.subplots()
            for outlier in outliers:
                plt.scatter(outlier.injection_order, outlier.tic, color="red")
                ax.annotate(outlier.id, (outlier.injection_order, outlier.tic))
            for spectrum in self._spectrum:
                if spectrum not in outliers:
                    plt.scatter(
                        spectrum.injection_order, spectrum.tic, color="blue")
                    ax.annotate(spectrum.id,
                                (spectrum.injection_order, spectrum.tic))
            for t in range(threshold):
                lbl = "%s times MAD" % (t + 1)
                ax.axhline(
                    mad * (t + 1), linestyle="--", color="rgboy" [t],
                    label=lbl
                    )
            plt.ylabel("Injection Order")
            plt.xlabel("Total Ion Count")
            plt.tight_layout()
            plt.legend(loc="best")
            plt.show()

        mad = _calculate_mad()
        outliers = _find_outliers(mad)

        if plot is True:
            _plot(outliers)

        if inplace is True:
            self._outlier_detected = True
            self.delete(outliers)
        else:
            return self.delete(outliers, inplace=False)

    def binning(self,
                bin_size=0.25,
                int_statistic="max",
                mass_statistic="mean",
                inplace=True):
        """Perform mass-binning.

        Parameters
        ----------

        bin_size : float, optional (default=0.25)
            Mass-to-ion bin size to use for binning.

        int_statistic : string, optional (default="median")
            The method used to calculate the binned intensitiy value.

            - If mean, calculated as the mean all spectrum intensity values for
              each bin.
            - If median, calculated as the median of all spectrum intensity
              values for each bin.
            - If count then compute the count of intensities within each bin.
            - If sum then compute the sum of intensity values within each bin.
            - If min then compute the minimum intensity value within each bin.
            - If max then compute the maximum intensity value within each bin.
        mass_statistic : string, optional (default="max")
            The method used to calculate the binned mass value.

            - If mean, calculated as the mean all spectrum mass values for
              a given bin.
            - If median, calculated as the median of all spectrum mass values
              for a given bin.

        inplace : boolean, optional (default=True)
            If False then return the binned Spectrums, else make the
            change within the object.

        """

        def _generate_ranged_bins():
            mass_range = self.get_mass_range()
            return np.arange(mass_range[0], mass_range[1], step=bin_size)

        def _calculate_mass_values(bins):
            mass_values = {b: [] for b in bins}
            for spectrum in self._spectrum:
                m = spectrum.masses
                for b in bins:
                    bindx = np.logical_and(m >= b, m <= (b + bin_size))
                    mass_values[b].extend(m[bindx])
            calculated_mass_values = []
            for bin, values in mass_values.iteritems():
                if values == []:
                    calculated_mass_values.append(bin)
                else:
                    method = getattr(np, mass_statistic)
                    calculated_mass_values.append(method(values))
            return sorted(mass_values)

        def _apply_binning(spectrum, bins):
            b_i, b_m, b_n = binned_statistic(
                spectrum.masses,
                spectrum.intensities,
                bins=bins,
                statistic=int_statistic)
            indx = np.invert(np.isnan(b_i))

            return b_m[:-1][indx], b_i[indx]

        bins = _generate_ranged_bins()

        if mass_statistic is not None:
            bins = _calculate_mass_values(bins)

        if inplace is False:
            t_sl = SpectrumList()

        for spectrum in self._spectrum:
            bins, intensities = _apply_binning(spectrum, bins)
            if inplace:
                spectrum.masses = bins
                spectrum.intensities = intensities
            else:
                t_s = copy(spectrum)
                t_s.masses = bins
                t_s.intensities = intensities
                t_sl.append(t_s)
        if inplace is False:
            return t_sl
        else:
            self._binned = True

    def normalise(self, method="tic", inplace=True):
        """Helper method to apply normalisation across all Spectrum objects
        within the SpectrumList.

        Normalisation aims to remove sources of variability within the spectrum

        Parameters
        ---------

        method : string, optional (default="tic")
            Method to use for normalisation.

            - If "tic" then apply TIC normalisation.
            - If "median" then apply normalisation by the median.

        inplace : boolean, optional (default=True)
            If False then return normalised intensities, else make the
            change within the object.

        """

        if inplace is True:
            for spectrum in self.tolist():
                spectrum.normalise(method=method)
            self._normalised = True
        else:
            t_sl = copy(self)
            for spectrum in t_sl.tolist():
                spectrum.normalise(method=method)
            return t_sl

    def transform(self, method="nlog", inplace=True):
        """Helper method to apply transformation across all Spectrum objects
        within the SpectrumList.

        Transformation aims to make the data less skewed.

        Parameters
        ----------

        method : string, optional (default="log10")
            Method to use for transformation.

            - If "log10" then apply log10 transformation.
            - If "cube" then apply cube transformation.
            - If "nglog" then apply nlog transformation.
            - If "log2" then apply log2 transformation.
            - If "glog" then apply globalised log transformation.
            - If "sqrt" then apply square root-based transformation.
            - If "ihs" then apply inverse hyperbolic sine transformation.

        inplace : boolean, optional (default=True)
            If False then return normalised intensities, else make the
            change within the object.

        """

        if inplace is True:
            for spectrum in self.tolist():
                spectrum.transform(method=method)
            self._transformed = True
        else:
            t_sl = copy(self)
            for spectrum in t_sl.tolist():
                spectrum.transform(method=method)
            return t_sl

    def scale(self, method="mc", inplace=True):
        """Apply mass scaling over the SpectrumList.

        Parameters
        ---------

        method : string, optional (default="mc")
            Method for applying mass-scaling.

            - If "mc" then apply mean centered scaling.
            - If "range" then apply ranged scaling.
            - If "auto" then apply auto scaling.
            - If "pareto" then apply pareto scaling.

        inplace : boolean, optional (default=True)
            If False then return a scaled SpectrumList, else make the
            change within the object.

        """

        def _mean_center(i):
            return i - np.mean(i)

        def _pareto(i):
            return _mean_center(i) / np.sqrt(np.std(i))

        def _range(i):
            return _mean_center(i) / (np.max(i) - np.min(i))

        def _auto(i):
            return _mean_center(i) / np.std(i)

        df = self.to_dataframe()

        for m in df:
            i = df[m].values
            if method.upper() == "MC":
                s_i = _mean_center(i)
            elif method.upper() == "RANGE":
                s_i = _range(i)
            elif method.upper() == "AUTO":
                s_i = _auto(i)
            elif method.upper() == "PARETO":
                s_i = _pareto(i)
            else:
                raise NotImplementedError(
                    "%s is not a valid scaler method" % method)
            df[m] = s_i

        sl = self.dataframe_to_spectrum_list(df)

        if inplace is True:
            self._scaled = True
            self._spectrum = sl._spectrum
        else:
            return sl

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
            raise ValueError(
                "%s is not a valid object type, need Spectrum object" %
                type(s))

    def value_imputation(self, method="basic", threshold=0.5, inplace=True):
        """Value imputator to replace missing values, and remove masses of which
        aren't imputatable.

        Parameters
        ---------

        method : string, optional (default="basic")
            Method for applying mass imputation.

            - If "basic" then imputate half the minimum intensity for the given
              spectrum.
            - If "mean" then imputate using the mean intensity for the given
              spectrum.
            - If "min" then imputate using the minimum intensity for the given
              spectrum.
            - If "median" then imputate using the median intensity for the
              given spectrum.

        threshold : float, optional (default=0.5)
            The threshold for the minimum number of missing intensities within
            a given mass.

        inplace : boolean, optional (default=True)
            If False then return a value imputated SpectrumList, else make the
            change within the object.

        """
        def _remove_by_threshold(df):
            null_count = df.isnull().sum()
            _t = len(df.index.values) * threshold
            to_keep = null_count <= threshold
            return df[df.columns[to_keep]]

        def _apply_imputation(df):
            for identifier in df.index:
                i = df.ix[identifier]
                if method.upper() == "BASIC":
                    filler = np.nanmin(i) / 2
                elif method.upper() == "MEAN":
                    filler = np.mean(i)
                elif method.upper() == "MIN":
                    filler = np.nanmin(i)
                elif method.upper() == "MEDIAN":
                    filler = np.nanmedian(i)
                else:
                    raise ValueError(
                        "%s is not a valid imputation method." % method
                        )
                df.ix[identifier] = df.ix[identifier].replace(np.nan, filler)
            return df

        df = self.to_dataframe()

        if method.upper() == "ALL":
            threshold = 0

        df = _remove_by_threshold(df)
        df = _apply_imputation(df)

        sl = self.dataframe_to_spectrum_list(df)
        if inplace is True:
            self._spectrum == sl._spectrum
            self._value_imputated = True
        else:
            return sl

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
        """Remove given Spectrum from the SpectrumList.

        Parameters
        ----------

        s : Spectrum
            A spectrum object or a list of spectrum objects to be removed.

        inplace : boolean, optional (default=True)
            If False then return the corrected intensities array, else make the
            change within the object.

        """
        if inplace is True:
            if type(s) == Spectrum:
                self._spectrum.remove(s)
            elif isinstance(s, list):
                [self._spectrum.remove(x) for x in s]

        else:
            cs = copy(self._spectrum)
            if type(s) == Spectrum:
                cs.remove(s)
            elif isinstance(s, list):
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

    def to_dataframe(self):
        """Flatten the SpectrumList to a Pandas DataFrame.

        Returns
        -------

        y : Pandas DataFrame
            A Pandas DataFrame containing the Spectrum Objects.

        """
        output = []
        for spectrum in self.tolist():
            df = pd.DataFrame(spectrum.intensities).T
            df.columns = spectrum.masses
            df.index = [spectrum.id]
            output.append(df)
        return pd.concat(output, axis=0)

    def tolist(self):
        """Return the SpectrumList as a list.

        Returns
        -------
        y : list
            A list containing Spectrum objects.

        """
        return self._spectrum

    def dataframe_to_spectrum_list(self, df):
        """Convert a Pandas DataFrame into a SpectrumList object.

        Parameters
        ----------

        df : Pandas DataFrame
            A Pandas DataFrame where the columns header are the masses,
            and each row denotes an individual spectrum.

        Returns
        -------
        df : Pandas DataFrame
            SpectrumList object containing the loaded Spectrum Objects.

        """

        t_sl = SpectrumList([])
        for id, values in df.iterrows():
            intensities = values.values
            idx = intensities != np.nan
            t_s = Spectrum(id=id)
            t_s.masses = df.columns[idx]
            t_s.intensities = intensities[idx]
            t_sl.append(t_s)
        return t_sl

    def get_mass_range(self):
        """Calculate the mass range of the entirety of the SpectrumList.

        Returns
        -------

        smallest : float
            The smallest mass value within the SpectrumList.

        largest : float
            The largest mass value within the SpectrumList.

        """
        smallest = None
        largest = None

        for spectrum in self.tolist():
            if min(spectrum.masses) < smallest or smallest is None:
                smallest = min(spectrum.masses)
            if max(spectrum.masses) > largest or largest is None:
                largest = max(spectrum.masses)

        return smallest, largest
