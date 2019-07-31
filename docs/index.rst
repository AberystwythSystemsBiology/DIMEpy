.. DIMEpy documentation master file, created by
   sphinx-quickstart on Thu Jul 18 11:57:51 2019.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to DIMEpy's documentation!
==================================

Python package for the high-throughput nontargeted metabolite
fingerprinting of nominal mass direct injection mass spectrometry
directly from mzML files.

This work is very much inspired by the methods detailed in
`High-throughput, nontargeted metabolite fingerprinting using nominal
mass flow injection electrospray mass spectrometry (Beckmann, et al,
2008) <https://www.nature.com/articles/nprot.2007.500>`__.

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

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   installation.rst
   getting_started.rst
   modules.rst

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
