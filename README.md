# DIMEpy: Direct Infusion MEtablomics (DIME) Processing in Python

**HERE BE DRAGONS:** This project is largely undocumented and untested, I do aim on sorting it all out eventually.

Python package for the high-thoroughput nontargeted metabolite fingerprinting of nominal mass direct injection mass spectrometry from ```mzML``` files.

Implementation of the methods detailed in:

```
High-throughput, nontargeted metabolite fingerprinting using nominal mass flow injection electrospray mass spectrometry

Beckmann, et al. (2008) - doi:10.1038/nprot.2007.500
```

## Installation

DIMEpy requires Python 2.7.+ and is unfortunately not compatible with Python 3.

You can install it through ```pypi``` using ```pip```:

```
pip install dimepy
```

Or alternatively install it manually using ```git```:

```
git clone https://www.github.com/KeironO/DIMEpy
cd DIMEpy
python setup.py install
```

## Bug reporting

Please report all bugs you find in the issues tracker. We would welcome all sorts of contribution, so please be as candid as you want.

## Contributors

* Keiron O'Shea (keo7@aber.ac.uk)

## Usage

The following script takes a path containing mzML files, processes them following the Beckmann, et al protocol and exports the result to an Excel file.


```python

# Importing modules required to run this script.

import dimepy
import os


# Path containing mzML files.
mzMLpaths = "/dir/to/mzMLs/"

# Path to save the output.
output_directory = "/output/directory/"


# mzML parameters.
parameters = {
    "MS1 Precision" : 1e-6,
    "MSn Precision" : 1e-6,
    "Measured Precision" : 1e-6,
    "Scan Range" : "apex",
    "Peak Type" : "peaks"
}


# Object to store processed spectrum.

for polarity in ["negative", "positive"]:
    spectrum_list = dimepy.SpectrumList()

    for index, file in enumerate(os.listdir(mzMLpaths)):
        # Read a mzML file from a given directory, and process it using given parameters.
        spectrum = dimepy.Spectrum(file_path=os.path.join(mzMLpaths, file),
                                   polarity=polarity, parameters=parameters)
        # Applying TIC normalisation
        spectrum.normalise(method="tic")
        # Applying generalised log transformation.
        spectrum.transform(method="glog")
        # Adding the processed spectrum to the spectrum list.
        spectrum_list.add(spectrum)


    # Create a spectrum list processor.
    processor = dimepy.SpectrumListProcessor(spectrum_list)

    # Apply MAD outlier detection.
    processor.outlier_detection()

    # Bin the spectrum to 0.25 m/z widths.
    processor.binning(bin_size=0.25)

    # Value imputation and value thresholding.
    processor.value_imputation(method="basic", threshold=0.5)

    # Applying mass-wise pareto scaling to the spectrum list.
    processor.scale(method="pareto")

    # Export the processed spectrum list back to a spectrum list object.
    processed_spectrum_list = processor.to_spectrumlist()

    # Convert the spectrum list to a Pandas DataFrame.
    df = processed_spectrum_list.flatten_to_dataframe()

    # Export processed spectrum to to Excel.
    df.to_excel(os.path.join(output_directory, polarity+".xlsx"))
```

## License

DIMEpy is licensed under the GNU General Public License v2.0.
