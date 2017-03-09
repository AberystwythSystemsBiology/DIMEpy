import numpy as np, pandas as pd, pickle as pkl, collections

class SpectrumList(object):
    def __init__(self, spectrum_list=[], processing_dict={}):
        self.spectrum_list = spectrum_list

    processor_dict = collections.OrderedDict()

    def __repr__(self):
        return self.spectrum_list

    def to_list(self):
        return list(self.spectrum_list)

    def append(self, spectrum):
        self.spectrum_list.append(spectrum)

    def remove(self, spectrum):
        self.spectrum_list.remove(spectrum)

    def pickle(self, fp="/tmp/pickled.pkl"):
        with open(fp, "wb") as output:
            pkl.dump(self, output, pkl.HIGHEST_PROTOCOL)

    def from_pickle(self, fp="/tmp/pickled.pkl"):
        with open(fp, "rb") as input:
            o = pkl.load(input)

        self.spectrum_list = o.spectrum_list
        self.processing_dict = o.processing_dict

    def to_csv(self, fp="/tmp/output.csv", delim=","):
        output = []
        for spectrum in self.spectrum_list:
            output.append([spectrum.id] + spectrum.masses.tolist())
            output.append([" "] + spectrum.intensities.tolist())
        pd.DataFrame(output).T.dropna().to_csv(fp, delimiter=delim, header=False, index=False)