import numpy as np, pandas as pd, pickle as pkl

class SpectrumList(object):
    def __init__(self, spectrum_list=[], binned=False):
        self.spectrum_list = spectrum_list
        self.binned = binned

    def __repr__(self):
        return repr(self.spectrum_list)

    def to_list(self):
        return list(self.spectrum_list)

    def add(self, spectrum):
        self.spectrum_list.append(spectrum)

    def remove(self, spectrum):
        self.spectrum_list.remove(spectrum)

    def outlier_detection(self, threshold=3):
        total_ion_counts = [np.nansum(s.intensities) for s in self.spectrum_list]
        mean_tic = np.nanmean(total_ion_counts)
        mean_ad= np.nanmean([abs(x-mean_tic) for x in total_ion_counts])

        ad_f_m = [abs((x-mean_tic)/mean_ad) for x in total_ion_counts]

        outlier_spectrums = [s for i, s in enumerate(self.spectrum_list) if ad_f_m[i] > threshold]

        for spectrum in outlier_spectrums:
            self.remove(spectrum)

    def bin(self, bw=0.25, removena=True):
        if self.binned == True:
            return
        for spectrum in self.spectrum_list:
            spectrum.bin(bw, removena)
            self.binned = True

    def smooth(self, sigma=1):
        for spectrum in self.spectrum_list:
            spectrum.smooth(sigma)

    def correct_baseline(self, lambda_=100, porder=1, max_iterations=15):
        for spectrum in self.spectrum_list:
            spectrum.correct_baseline(lambda_, porder, max_iterations)

    def normalise(self, method="tic"):
        for spectrum in self.spectrum_list:
            spectrum.normalise(method)

    def transform(self, method="log10"):
        for spectrum in self.spectrum_list:
            spectrum.transform(method)

    def pickle(self, fp="/tmp/pickled.pkl"):
        with open(fp, "wb") as output:
            pkl.dump(self, output, pkl.HIGHEST_PROTOCOL)

    def from_pickle(self, fp="/tmp/pickled.pkl"):
        with open(fp, "rb") as input:
            o = pkl.load(input)

        self.spectrum_list = o.spectrum_list
        self.binned = o.binned

    def to_excel(self, fp="/tmp/output.xlsx"):
        if self.binned == True:
            output = []
            for spectrum in self.spectrum_list:
                df = pd.DataFrame(spectrum.intensities).T
                df.columns = spectrum.masses
                df.index = [str(spectrum.identifier)]
                output.append(df)
            df = pd.concat(output, axis=0)
            df.to_excel(fp)
        else:
            output = []
            for spectrum in self.spectrum_list:
                output.append([spectrum.identifier] + spectrum.masses.tolist())
                output.append([""] + spectrum.intensities.tolist())
            pd.DataFrame(output).T.dropna().to_excel(fp, index=False, header=False)

    def to_csv(self, fp="/tmp/output.csv", delim=","):
        output = []
        for spectrum in self.spectrum_list:
            output.append([spectrum.identifier] + spectrum.masses.tolist())
            output.append([" "] + spectrum.intensities.tolist())
        pd.DataFrame(output).T.dropna().to_csv(fp, delimiter=delim, header=False, index=False)