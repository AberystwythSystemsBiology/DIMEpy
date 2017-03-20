import pandas as pd, pickle as pkl, collections, csv, numpy as np

class SpectrumList(object):
    def __init__(self, spectrum_list=[], processing_dict=collections.OrderedDict()):
        self.spectrum_list = spectrum_list
        self.processing_dict = processing_dict

    def __repr__(self):
        return self.spectrum_list

    def to_list(self):
        return list(self.spectrum_list)

    def append(self, spectrum):
        if spectrum.id not in [x.id for x in self.spectrum_list]:
            self.spectrum_list.append(spectrum)
        else:
            # TODO: THROW AN EXCEPTION
            return

    def remove(self, spectrum):
        self.spectrum_list.remove(spectrum)

    def pickle(self, fp="/tmp/pickled.pkl"):
        with open(fp, "wb") as output:
            pkl.dump(self, output, pkl.HIGHEST_PROTOCOL)

    def from_pickle(self, fp="/tmp/pickled.pkl"):
        with open(fp, "rb") as input:
            o = pkl.load(input)

        self.spectrum_list = o.spectrum_list
        try:
            # Edit this to fit the actual init
            self.processing_dict = o.processor_dict
        except:
            self.processing_dict = None

    def to_csv(self, fp="/tmp/output.csv", delim=","):
        output = []
        for spectrum in self.spectrum_list:
            output.append([spectrum.id] + spectrum.masses.tolist())
            output.append([" "] + spectrum.intensities.tolist())
        output = pd.DataFrame(output).T.replace(np.nan, "", regex=True).values
        with open(fp, "wb") as out_file:
            writer = csv.writer(out_file, delimiter=delim)
            for row in output:
                writer.writerow(row)
            out_file.close()

    def from_csv(self):
        # TODO
        pass

    def to_excel(self, fp="/tmp/output.xlsx"):
        if "binning" or "center" in self.processor_dict.keys():
            output = []
            for spectrum in self.spectrum_list:
                df = pd.DataFrame(spectrum.intensities).T
                df.columns = spectrum.masses
                df.index = [spectrum.id]
                output.append(df)
            output = pd.concat(output, axis=0)
            output.to_excel(fp)
        else:
            return

