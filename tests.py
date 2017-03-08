from dimspy.Spectrum import Spectrum
from dimspy.SpectrumList import SpectrumList

# File directory where all of the mzML files are found.
file_directory = "/home/keo7/Data/mzMLs/"


spectrum_list = SpectrumList()
spectrum_list.from_pickle("test.pkl")


spectrum_list.smooth()
spectrum_list.correct_baseline()
spectrum_list.bin(0.25)
spectrum_list.value_imputation()
spectrum_list.transform()
spectrum_list.normalise()

import matplotlib.pyplot as plt

plt.figure()
for spectrum in spectrum_list.to_list():
    plt.plot(spectrum.masses, spectrum.intensities, label=spectrum.identifier)
plt.legend()
plt.show()

'''
# Use os to list all the files in the file directory
import os
file_list = os.listdir(file_directory)

# Create a new SpectrumList to store Spectrum objects
spectrum_list = SpectrumList()

# Loop through the file list, you could probably make it faster with joblib.
for index, file_name in enumerate(file_list):
    print index+1, "/", len(file_list)
    # Using the file_name as an identifier.
    sample_id = file_name.split(".")[0]
    # Creating a blank Spectrum object, ensuring that the identifier isn't blank.
    spectrum = Spectrum(identifier=sample_id)
    # Filling the spectrum object with data from a mzML file.
    spectrum.from_mzml(filepath=file_directory+file_name,
                                  polarity="positive", scan_range="apex", peak_type="peaks")

    # Add the Spectrum object to the SpectrumList
    spectrum_list.add(spectrum)


spectrum_list.pickle("test.pkl")

'''