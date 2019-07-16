# DIMEpy: Direct Infusion MEtablomics processing in python

Python package for the high-thoroughput nontargeted metabolite fingerprinting of nominal mass direct injection mass spectrometry directly from ```mzML``` files.

This is very much an implementation and extension of the  of the methods detailed in High-throughput, nontargeted metabolite fingerprinting using nominal mass flow injection electrospray mass spectrometry (Beckmann, et al, 2008).

## Features

- Loading mass spectrometry files from mzML.
  - Support for polarity switching.
  - MAD-estimated infusion profiling.
- Spurious peak elimination.
- Spectrum export for direct dissemination using Metaboanalyst.
- Spectral binning.
- Spectral normalisation.
  - TIC
  - Median
  - Mean
- Spectral transformation.
  - log10
  - cube
  - nlog
  - log2
  - glog
  - sqrt
  - ihs

## Usage

```python
# This will be rewritten shortly.
```

## Installation

DIMEpy requires Python 3+ and is unfortunately not compatible with Python 2. If you are still using Python 2, a clever workaround is to install Python 3 and use that instead.

You can install it through ```pypi``` using ```pip```:

```
pip install dimepy
```

If you want the 'bleeding edge' version this, you can also install directly from this repository using ```git```:

```
git clone https://www.github.com/AberystwythSystemsBiology/DIMEpy
cd DIMEpy
python setup.py install
```

Or use ```git``` and ```pip``` in unison:

```
pip install git+https://www.github.com/AberystwythSystemsBiology/DIMEpy
```

## Bug reporting

Please report all bugs you find in the issues tracker, please do not email me directly as I'm struggling to keep track of what needs to be fixed. We welcome all sorts of contribution, so please be as candid as you want.

## Contributors

* **Lead Developer:** Keiron O'Shea (keo7@aber.ac.uk)
* **Project Supervisor:** Chuan Lu (cul@aber.ac.uk)
* **Project Supervisor:** Luis AJ Mur (lum@aber.ac.uk)
* **Methods Expert:** Manfred Beckmann (meb@aber.ac.uk)

## License

DIMEpy is licensed under the GNU General Public License v2.0.
