.. DIMEpy documentation master file, created by
   sphinx-quickstart on Thu Jul 18 11:57:51 2019.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to DIMEpy's documentation!
==================================

|PyPI - Python Version| |PyPI| |PyPI - License| |DOI| |PyPI - Status|

Python package for the high-throughput nontargeted metabolite
fingerprinting of nominal mass direct injection mass spectrometry
directly from mzML files.

Features
--------

-  Loading mass spectrometry files from mzML.

   -  Support for polarity switching.
   -  MAD-estimated infusion profiling.

-  Assay-wide outlier spectrum detection.
-  Spurious peak elimination.
-  Spectrum export for direct dissemination using Metaboanalyst.
-  Spectral binning.
-  Value imputation.
-  Spectral normalisation.

   -  including TIC, median, mean…

-  Spectral transformation.

   -  including log10, cube, nlog, log2, glog, sqrt, ihs…

-  Export to array for statistical analysis in Metaboanalyst.

Contributors
------------

-  **Lead Developer:** Keiron O’Shea (keo7@aber.ac.uk)
-  **Developer:** Rob Bolton (rab26@aber.ac.uk)
-  **Project Supervisor:** Chuan Lu (cul@aber.ac.uk)
-  **Project Supervisor:** Luis AJ Mur (lum@aber.ac.uk)
-  **Methods Expert:** Manfred Beckmann (meb@aber.ac.uk)

License
-------

DIMEpy is licensed under the `GNU General Public License
v3.0 <https://raw.githubusercontent.com/AberystwythSystemsBiology/DIMEpy/master/LICENSE>`__.

.. |PyPI - Python Version| image:: https://img.shields.io/pypi/pyversions/DIMEpy.svg
.. |PyPI| image:: https://img.shields.io/pypi/v/DIMEpy.svg
.. |PyPI - License| image:: https://img.shields.io/pypi/l/DIMEpy.svg
.. |DOI| image:: https://zenodo.org/badge/DOI/10.5281/zenodo.3340120.svg
.. |PyPI - Status| image:: https://img.shields.io/pypi/status/DIMEpy.svg

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   installation.rst
   getting_started.rst
   modules.rst
   example_scripts.rst

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
