from dimspy.Spectrum import Spectrum
from dimspy.SpectrumList import SpectrumList


# Creating a new spectrum object with custom masses and intensities
k = Spectrum(identifier="A10", masses=[10, 20, 30], intensities=[3.43, 3.54, 3.21])

# File directory where all of the mzML files are found.
file_directory = "/home/keo7/Data/mzMLs/"

# Use os to list all the files in the file directory
import os
file_list = os.listdir(file_directory)

# Create a new SpectrumList to store Spectrum objects
spectrum_list = SpectrumList()
spectrum_list.from_pickle("out.pkl")

# Loop through the file list, you could probably make it faster with joblib.
for index, file_name in enumerate(file_list):
    if index > 0:
        break
    # Using the file_name as an identifier.
    sample_id = file_name.split(".")[0]
    # Creating a blank Spectrum object, ensuring that the identifier isn't blank.
    spectrum = Spectrum(identifier=sample_id)
    # Filling the spectrum object with data from a mzML file.
    spectrum.from_mzml(filepath=file_directory+file_name,
                                  polarity="positive", scan_range="apex", peak_type="peaks")


    print spectrum.identifier
    # Add the Spectrum object to the SpectrumList
    spectrum_list.add(spectrum)



'''
spectrum_list.from_pickle("out.pkl")

print spectrum_list

# Saving the raw non-processed values to a pickle file.
spectrum_list.pickle("out.pkl")
'''

import matplotlib.pyplot as plt

plt.plot()
for spectrum in spectrum_list.to_list():
    plt.plot(spectrum.masses, spectrum.intensities, label=spectrum.identifier)
plt.legend()
plt.show()