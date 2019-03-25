# DIMEpy: Direct Infusion MEtablomics (DIME) Processing in Python

Python package for the high-thoroughput nontargeted metabolite fingerprinting of nominal mass direct injection mass spectrometry from ```mzML``` files.

Implementation of the methods detailed in:

```
High-throughput, nontargeted metabolite fingerprinting using nominal mass flow injection electrospray mass spectrometry

Beckmann, et al. (2008) - doi:10.1038/nprot.2007.500
```

## Installation

DIMEpy requires Python 3+ and is unfortunately not compatible with Python 2. If you are still using Python 2, a clever workaround is to install Python 3 and use that instead.

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

The following script takes a path containing mzML files, processes them following the Beckmann, et al protocol and exports the result to a comma-seperated file.

```python

# Importing modules required to run this script.
from dimepy import Spectrum, SpectrumList
import os

# Path containing mzML files.
mzMLpaths = "/dir/to/mzMLs/"

# List to store the spectra
l = []

for index, file in enumerate(os.listdir(mzMLpaths)):
  filepath = os.path.join(mzMLpaths, file)
  # Load in the spectrum directly using default parameters.
  spectrum = Spectrum(filepath, identifier=file)
  spec.load()
  spec.get_polarity("positive")
  spec.get_apex(3)
  spec.get()

  for scan in spec.scans:
    scan.bin()

  spec.remove_spurious_peaks()
  spec.get()

  l.append(spec)

spectrum_list = SpectrumList(l, bin_width=0.02)

spectrum_list.df.to_csv("/tmp/raw_ish.csv")

spectrum_list.transform()
spectrum_list.normalise()

spectrum_list.value_imputate()

spectrum_list.df.to_csv("/tmp/final.csv")

```

## License

DIMEpy is licensed under the GNU General Public License v2.0.
