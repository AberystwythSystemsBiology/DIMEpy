'''
    DIMSpy - example.py


    An example Python script, demonstrating some of the functionality of the package.

    Author: Keiron O'Shea <keo7@aber.ac.uk>

'''

from dimspy.Spectrum import Spectrum
from dimspy.SpectrumList import SpectrumList

# Creating an example Spectrum object
s1 = Spectrum(identifier="Example", masses=[1, 2, 3, 4], intensities=[4, 3, 2, 1])

# Printing out the Spectrum ID.
print s1.identifier

s1.smooth(sigma=4)
print "Smoothed intensities:", s1.intensities

# Creating a Spectrum List object, for multi-sample processing
sl = SpectrumList()

# Adding our Spectrum object to the SpectrumList
sl.add(s1)

# Returning Spectrum objects as a list
print sl.to_list()



# Smoothing over all Spectrum objects in the SpectrumList
sl.smooth()
'''
    Saving our list in a human readable Excel file.

    The output differs on the type of SpectrumList is being passed.

    Non-binned

    |Identifier|         |
    |   mass   |intensity|

    Binned
    |         | mass1 | mass2 |....
    |Idenifier| int1  |  int2 |....
'''
sl.to_excel("/tmp/out.xlsx")


# Pickling our object for later use.
sl.pickle("/tmp/pickled_list.pkl")


# Creating a blank SpectrumList
s2 = SpectrumList()
# Loading our pickled SpectrumList
s2.from_pickle("/tmp/pickled_list.pkl")