# DIMEpy: Direct Infusion MEtablomics processing in python

[![Project Status: Active - The project has reached a stable, usable state and is being actively developed.](http://www.repostatus.org/badges/0.1.0/active.svg)](http://www.repostatus.org/#active)
![PyPI - Python Version](https://img.shields.io/pypi/pyversions/DIMEpy.svg)
![PyPI](https://img.shields.io/pypi/v/DIMEpy.svg)
![PyPI - License](https://img.shields.io/pypi/l/DIMEpy.svg)
![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3340120.svg)
![PyPI - Status](https://img.shields.io/pypi/status/DIMEpy.svg)
![GitHub commit activity](https://img.shields.io/github/commit-activity/y/KeironO/DIMEpy.svg)

Python package for the high-thoroughput nontargeted metabolite fingerprinting of nominal mass direct injection mass spectrometry directly from mzML files.

This work is very much inspired by the methods detailed in [High-throughput, nontargeted metabolite fingerprinting using nominal mass flow injection electrospray mass spectrometry (Beckmann, et al, 2008)](https://www.nature.com/articles/nprot.2007.500). 

## Features

- Loading mass spectrometry files from mzML.
  - Support for polarity switching.
  - MAD-estimated infusion profiling.
- Assay-wide outlier spectrum detection.
- Spurious peak elimination.
- Spectrum export for direct dissemination using Metaboanalyst.
- Spectral binning.
- Value imputation.
- Spectral normalisation.
  - including TIC, median, mean...
- Spectral transformation.
  - including log10, cube, nlog, log2, glog, sqrt, ihs...
- Export to array for statistical analysis in Metaboanalyst.

## Installation

DIMEpy requires Python 3+ and is unfortunately not compatible with Python 2. If you are still using Python 2, a clever workaround is to install Python 3 and use that instead.

You can install it through ```pypi``` using ```pip```:

```
pip install dimepy
```

If you want the 'bleeding edge' version this, you can also install directly from this repository using ```git``` - but beware of dragons:

```
pip install git+https://www.github.com/AberystwythSystemsBiology/DIMEpy
```

## Usage

To use the package, type the following into your Python console:

```python
>>> import dimepy
```

At the moment, this pipeline only supports mzML files. You can easily convert proprietary formats to mzML using [ProteoWizard](http://www.proteowizard.org/download.html).

### Loading a single file

If you are only going to load in a single file for fingerprint matrix estimation, then just create a new spectrum object. If the sample belongs to a characteristic, it is recommend that you also pass it through when instantiating a new ```Spectrum``` object.

```python
>>> filepath = "/file/to/file.mzML"
>>> spec = dimepy.Spectrum(filepath, identifier="example", stratification="class_one")
/file/to/file.mzML
```

By default the Spectrum object doesn't set a snr estimator. It is **strongly recommended** that you set a signal to noise estimation method when instantiating the Spectrum object.

If your experimental protocol makes use of mixed-polarity scanning, then please ensure that you limit the scan ranges to best match what polarity you're interested in analysing:

```python
>>> spec.limit_polarity("negative")
```


If you are using FIE-MS it is strongly recommended that you use just the infusion profile to generate your mass spectrum. For example, if your scan profiles look like this:

```
        |        _
      T |       / \
      I |      /   \_
      C |_____/       \_________________
        0     0.5     1     1.5     2 [min]
```

Then it is fair to assume that the infusion occured during the scans ranging from 30 seconds to 1 minute. The ```limit_infusion()``` method does this by estimating the mean absolute deviation (MAD) of total ion counts (TIC) before limiting the profile to the range between the time range in which whatever multiple of MAD has been estimated:

```python
>>> spec.limit_infusion(2) # 2 times the MAD.
```

Now, we are free to load in the scans to generate a base mass_spectrum:

```python
>>> spec.load_scans()
```

You should now be able to access the generated mass spectrum using the ```masses``` and ```intensities``` attributes:

```python
>>> spec.masses
array([ ... ])
>>> spec.intensities
array([ ... ])
```

### Working with multiple files

A more realistic pipeline would be to use multiple mass-spectrum files. This is where things really start to get interesting. The ```SpectrumList``` object facilitates this through the use of the ```append``` method:

```python
>>> speclist = dimepy.SpectrumList()
>>> speclist.append(spec)
```

You can make use of an iterator to recursively generate ```Spectrum``` objects, or do it manually if you want.

If you're only using this pipeline to extract mass spectrum for Metabolanalyst, then you can now simply call the ```_to_csv``` method:

```python
>>> speclist.to_csv("/path/to/output.csv", output_type="metaboanalyst")
```

That being said, this pipeline contains many of the preprocessing methods found in Metaboanalyst - so it may be easier for you to just use ours.

As a diagnostic measure, the TIC can provide an estimation of factos that may adversely affect the overal intensity count of a run. As a rule, it is common to remove spectrum in which the TIC deviates 2/3 times from the median-absolute deviation. We can do this by calling the ```detect_outliers``` method:

```python
>>> speclist.detect_outliers(thresh = 2, verbose=True)
Detected Outliers: outlier_one;outlier_two 
```

A common first step in the analysis of mass-spectrometry data is to bin the data to a given mass-to-ion value. To do this for all ```Spectrum``` held within our ```SpectrumList``` object, simply apply the ```bin``` method:

```python
>>> speclist.bin(0.25) # binning our data to a bin width of 0.25 m/z
```

In FIE-MS null values should concern no more than 3% of the total number of identified bins. However, imputation is required to streamline the analysis process (as most multivariate techniques are unable to accomodate missing data points). To perform value imputation, just use ```value_imputate```:

```python
>>> speclist.value_imputate()
```

Now transforming and normalisating the the spectrum objects in an samples independent fashion can be done using the following:

```python
>>> speclist.transform()
>>> speclist.normalise()
```

Once completed, you are now free to export the data to a data matrix:

```python
>>> speclist.to_csv("/path/to/proc_metabo.csv", output_type="matrix")
```

This should give you something akin to:

| Sample ID | M0 | M1 | M2 | M3 |... |
|-----------|----|----|----|----|----|
| Sample 1 | 213 | 634 | 3213 | 546 | ... |
| Sample 2 | 132 | 34 | 713 | 6546 |... |
| Sample 3 | 1337  | 42 | 69 | 420 | ... |

## Bug reporting and feature suggestions

Please report all bugs or feature suggestions to the [issues tracker](https://github.com/AberystwythSystemsBiology/DIMEpy/issues). **Please do not email me directly** as I'm struggling to keep track of what needs to be fixed. 

We welcome all sorts of contribution, so please be as candid as you want(!)

## Contributors

* **Lead Developer:** Keiron O'Shea (keo7@aber.ac.uk)
* **Project Supervisor:** Chuan Lu (cul@aber.ac.uk)
* **Project Supervisor:** Luis AJ Mur (lum@aber.ac.uk)
* **Methods Expert:** Manfred Beckmann (meb@aber.ac.uk)

## License

DIMEpy is licensed under the [GNU General Public License v2.0](https://raw.githubusercontent.com/AberystwythSystemsBiology/DIMEpy/master/LICENSE).
