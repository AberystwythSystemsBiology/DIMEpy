import os
import dimspy
import cPickle as pickle

data_directory = "/home/keo7/Data/Denisa_Saliva/mzMLs/"

# mzML parameters.
parameters = {
    "MS1 Precision" : 1e-6,
    "MSn Precision" : 1e-6,
    "Measured Precision" : 1e-6,
    "Scan Range" : "apex",
    "Peak Type" : "peaks"
}


spec_list = []


if os.path.exists("test.p") != True:
    for index, file in enumerate(os.listdir(data_directory)):
        if index > 20:
            break

        spectrum = dimspy.Spectrum(file_path=os.path.join(data_directory, file),
                                   polarity="negative", parameters=parameters)

        spectrum.transform()
        spectrum.normalise()

        spec_list.append(spectrum)
    processor = dimspy.SpectrumListProcessor(dimspy.SpectrumList(spec_list))
    pickle.dump(processor, open("test.p", "wb"))

else:
    processor = pickle.load(open("test.p", "rb"))

processor.outlier_detection()
processor.binning(bin_size=0.25)
processor.center()
processor.value_imputation()

sl = processor.to_spectrumlist()
print sl.flatten_to_dataframe()



