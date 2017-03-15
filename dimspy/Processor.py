import numpy as np, scipy.stats as stats, pandas as pd, collections, json

class SpectrumListProcessor(object):

    def __init__(self, spectrum_list):
        self.spectrum_list = spectrum_list

    processor_dict = collections.OrderedDict()
    binned_intensities = None

    def outlier_detection(self, threshold=3):
        total_ion_counts = [np.nansum(s.intensities) for s in self.spectrum_list.to_list()]
        mean_tic = np.nanmean(total_ion_counts)
        mean_ad = np.nanmean([abs(x - mean_tic) for x in total_ion_counts])

        ad_f_m = [abs((x - mean_tic) / mean_ad) for x in total_ion_counts]

        outlier_spectrums = [s for i, s in enumerate(self.spectrum_list.to_list()) if ad_f_m[i] > threshold]

        for spectrum in outlier_spectrums:
            self.spectrum_list.remove(spectrum)

        self.processor_dict["outlier_detection"] = {
            "removed" : [x.id for x in outlier_spectrums],
            "mean absolute deviation" : mean_ad,
            "threshold" : str(threshold)+" +/- MAD"
        }

    def smooth(self, sigma=3):
        for spectrum in self.spectrum_list.to_list():
            spectrum.smooth(sigma)

        self.processor_dict["smoothing"] = {
            "sigma" : sigma
        }

    def correct_baseline(self, lambda_=100, porder=1, max_iterations=15):
        for spectrum in self.spectrum_list.to_list():
            spectrum.correct_baseline(lambda_=100, porder=1, max_iterations=15)

        self.processor_dict["baseline correction"] = {
            "maximum iterations" : max_iterations
        }

    def peak_detection(self, delta=3):
        for spectrum in self.spectrum_list.to_list():
            spectrum.get_peaks(delta)

        self.processor_dict["peak detection"] = {
            "delta" : delta
        }

    def normalise(self, method="tic"):
        for spectrum in self.spectrum_list.to_list():
            spectrum.normalise(method)

        self.processor_dict["normalise"] = {
            "method" : method
        }

    def transform(self, method="log10"):
        for spectrum in self.spectrum_list.to_list():
            spectrum.transform(method)

        self.processor_dict["transformation"] = {
            "method" : method
        }

    def binning(self, bin_size=0.05):
        shared_masses = []

        def _bin(spectrum):
            bins = np.arange(round(min(spectrum.masses)), round(max(spectrum.masses)), step=bin_size)
            b_intensities, b_masses, b_num = stats.binned_statistic(spectrum.masses, spectrum.intensities, bins=bins)
            spectrum.masses = b_masses[:-1]
            spectrum.intensities = b_intensities
            return b_masses

        for spectrum in self.spectrum_list.to_list():
            binned_masses = _bin(spectrum)
            shared_masses.append(binned_masses.tolist())

        self.binned_intensities = np.array(sum(shared_masses, []))
        self.processor_dict["binning"] = {
            "bin size" : bin_size,
            "total number of bins" : self.binned_intensities.size
        }

    def value_imputation(self, method="knn", threshold=0.5):

        def _remove_bins_by_threshold():
            sample_threshold = len(self.spectrum_list.to_list()) * threshold
            output = []
            for spectrum in self.spectrum_list.to_list():
                df = pd.DataFrame(spectrum.intensities).T
                df.columns = spectrum.masses
                df.index = [spectrum.id]
                output.append(df)
            df = pd.concat(output, axis=0)
            df.dropna(axis=1, thresh=sample_threshold, inplace=True)
            return df

        def _value_imputation(df):
            if method == "knn":
                from sklearn.preprocessing import Imputer
                imp = Imputer(axis=0)
                imputated = imp.fit_transform(df)
                df = pd.DataFrame(imputated, columns=df.columns, index=df.index)
            elif method == "basic":
                df.fillna(value=(np.nanargmin(df.values) / 2), inplace=True)
            return df

        def _back_to_spectrum(df):
            for id, values in df.iterrows():
                masses = np.array(list(values.index))
                intensities = values.values
                spectrum = [x for x in self.spectrum_list.to_list() if x.id == id][0]
                spectrum.masses = masses
                spectrum.intensities = intensities


        df = _remove_bins_by_threshold()
        df = _value_imputation(df)

        _back_to_spectrum(df)

        self.binned_intensities = np.array(df.columns)

        self.processor_dict["value imputation"] = {
            "method" : method,
            "sample_threshold" : threshold
        }


    def to_spectrumlist(self):
        from SpectrumList import SpectrumList
        return SpectrumList(spectrum_list=self.spectrum_list.to_list(), processing_dict=self.processor_dict)

    def save_processor_dict(self, fp="/tmp/processor_dict.json"):
        with open(fp, "wb") as out_file:
            json.dump(self.processor_dict, out_file, indent=4)