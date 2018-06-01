# DIMEpy: Direct Infusion MEtablomics (DIME) Processing in Python

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

alternatively install it manually using ```git```:

```
git clone https://www.github.com/KeironO/DIMEpy
cd DIMEpy
python setup.py install
```

Or use ```git``` and ```pip``` in unison.

```
pip install git+https://www.github.com/KeironO/DIMEpy
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

# Where we'll store the spectrum.
spectrum_list = dimepy.SpectrumList()

for index, file in enumerate(os.listdir(mzMLpaths)):
  # Load in the spectrum directly using default parameters.
  spectrum = dimepy.Spectrum(os.path.join(mzMLpaths, file))
  # Correct for baseline.
  spectrum.baseline_correction(qtl=0.6)
  spectrum_list.append(spectrum)

# Write the raw spectrum to a comma seperated file.
spectrum_list.to_csv("raw.csv")
# Convert the object to a SpectrumListProcessor for processing.

# Apply outlier detection to remove spurious samples.
spectrum_list.outlier_detection()
# Bin masses over 0.125 m/z.
spectrum_list.binning(bin_size=0.125)
# Value imputate where < 50% of the values are lost across all samples.
spectrum_list.value_imputation(method="basic", threshold=0.5)
# Normalise over the total ion count.
spectrum_list.normalise(method="TIC")
# Apply generalised log transformation
spectrum_list.transform(method="glog")

# Write the processed spectrum to a comma seperated file.
spectrum_list.to_csv("processed.csv")
```

## License

DIMEpy is licensed under the GNU General Public License v2.0.
