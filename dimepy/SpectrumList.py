# -*- coding: utf-8 -*-

import pandas as pd
import cPickle as pkl
import collections
import csv
import numpy as np


class SpectrumList(object):
    def __init__(self, _spectrum=[]):
        self.__spectrum = _spectrum

    def append(self, Spectrum):
        self.__spectrum.append(Spectrum)

    def pickle(self, file_path):
        with open(file_path, "wb") as output:
            pkl.dump(self, output, pkl.HIGHEST_PROTOCOL)

    def remove(self, Spectrum):
        self.__spectrum.remove(Spectrum)

    def to_csv(self, file_path, delimiter=","):
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
        output = []
        for spectrum in self.to_list():
            df = pd.DataFrame(spectrum.intensities).T
            df.columns = spectrum.masses
            df.index = [spectrum.id]
            output.append(df)
        return pd.concat(output, axis=0)

    def to_list(self):
        return list(self.__spectrum)

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

        for spectrum in self.to_list():
            if min(spectrum.masses) < smallest or smallest is None:
                smallest = min(spectrum.masses)
            if max(spectrum.masses) > largest or largest is None:
                largest = max(spectrum.masses)

        return smallest, largest

    def __repr__(self):
        return repr(self.to_list())
