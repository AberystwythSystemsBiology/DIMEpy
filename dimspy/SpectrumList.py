import pandas as pd
import cPickle as pkl
import collections
import csv
import numpy as np


class SpectrumList(object):
    def __init__(self, _spectrum=[]):
        '''

        :param _spectrum:
        '''
        self.__spectrum = _spectrum

    def append(self, Spectrum):
        '''

        :param Spectrum:
        :return:
        '''
        self.__spectrum.append(Spectrum)

    def pickle(self, file_path):
        '''

        :param file_path:
        :return:
        '''
        with open(file_path, "wb") as output:
            pkl.dump(self, output, pkl.HIGHEST_PROTOCOL)

    def remove(self, Spectrum):
        '''

        :param Spectrum:
        :return:
        '''
        self.__spectrum.remove(Spectrum)

    def to_csv(self, file_path, delimiter=","):
        '''

        :param file_path:
        :param delimiter:
        :return:
        '''
        raw_values = []
        for spectrum in self.__spectrum:
            raw_values.append([spectrum.id] + spectrum.masses.tolist())
            raw_values.append([" "] + spectrum.intensities.tolist())
        raw_values = pd.DataFrame(raw_values).T.replace(
            np.nan, "", regex=True).values

        with open(file_path, "wb") as output:
            writer = csv.writer(output, delimiter=delimiter)
            for row in raw_values:
                writer.writerow(row)
            output.close()

    def flatten_to_dataframe(self):
        '''

        :return:
        '''
        output = []
        for spectrum in self.to_list():
            df = pd.DataFrame(spectrum.intensities).T
            df.columns = spectrum.masses
            df.index = [spectrum.id]
            output.append(df)
        return pd.concat(output, axis=0)

    def to_list(self):
        return list(self.__spectrum)

    def get_mass_range(self):
        smallest = None
        largest = None

        for spectrum in self.to_list():
            if min(spectrum.masses) < smallest or smallest is None:
                smallest = min(spectrum.masses)
            if max(spectrum.masses) > largest or largest is None:
                largest = max(spectrum.masses)

        return smallest, largest

    def __repr__(self):
        return repr(self.to_list())
