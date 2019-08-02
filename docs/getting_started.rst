Getting Started
===============

To use the package, type the following into your Python console:

.. code:: python

   >>> import dimepy

At the moment, this pipeline only supports mzML files. You can easily
convert proprietary formats to mzML using
`ProteoWizard <http://www.proteowizard.org/download.html>`__.

Loading a single file
~~~~~~~~~~~~~~~~~~~~~

If you are only going to load in a single file for fingerprint matrix
estimation, then just create a new spectrum object. If the sample
belongs to a characteristic, it is recommend that you also pass it
through when instantiating a new ``Spectrum`` object.

.. code:: python

   >>> filepath = "/file/to/file.mzML"
   >>> spec = dimepy.Spectrum(filepath, identifier="example", stratification="class_one")
   /file/to/file.mzML

By default the Spectrum object doesn’t set a snr estimator. It is
**strongly recommended** that you set a signal to noise estimation
method when instantiating the Spectrum object.

If your experimental protocol makes use of mixed-polarity scanning, then
please ensure that you limit the scan ranges to best match what polarity
you’re interested in analysing:

.. code:: python

   >>> spec.limit_polarity("negative")

If you are using FIE-MS it is strongly recommended that you use just the
infusion profile to generate your mass spectrum. For example, if your
scan profiles look like this:

::

           |        _
         T |       / \
         I |      /   \_
         C |_____/       \_________________
           0     0.5     1     1.5     2 [min]

Then it is fair to assume that the infusion occured during the scans
ranging from 30 seconds to 1 minute. The ``limit_infusion()`` method
does this by estimating the median absolute deviation (MAD) of total ion
count (TIC) before limiting the profile to the range between the time
range in which whatever multiple of MAD has been estimated:

.. code:: python

   >>> spec.limit_infusion(2) # 2 times the MAD.

Now, we are free to load in the scans to generate a base mass_spectrum:

.. code:: python

   >>> spec.load_scans()

You should now be able to access the generated mass spectrum using the
``masses`` and ``intensities`` attributes:

.. code:: python

   >>> spec.masses
   array([ ... ])
   >>> spec.intensities
   array([ ... ])

Working with multiple files
~~~~~~~~~~~~~~~~~~~~~~~~~~~

A more realistic pipeline would be to use multiple mass-spectrum files.
This is where things really start to get interesting. The
``SpectrumList`` object facilitates this through the use of the
``append`` method:

.. code:: python

   >>> speclist = dimepy.SpectrumList()
   >>> speclist.append(spec)

You can make use of an iterator to recursively generate ``Spectrum``
objects, or do it manually if you want.

If you’re only using this pipeline to extract mass spectrum for
Metabolanalyst, then you can now simply call the ``_to_csv`` method:

.. code:: python

   >>> speclist.to_csv("/path/to/output.csv", output_type="metaboanalyst")

That being said, this pipeline contains many of the preprocessing
methods found in Metaboanalyst - so it may be easier for you to just use
ours.

As a diagnostic measure, the TIC can provide an estimation of factos
that may adversely affect the overal intensity count of a run. As a
rule, it is common to remove spectrum in which the TIC deviates 2/3
times from the median-absolute deviation. We can do this by calling the
``detect_outliers`` method:

.. code:: python

   >>> speclist.detect_outliers(thresh = 2, verbose=True)
   Detected Outliers: outlier_one;outlier_two

A common first step in the analysis of mass-spectrometry data is to bin
the data to a given mass-to-ion value. To do this for all ``Spectrum``
held within our ``SpectrumList`` object, simply apply the ``bin``
method:

.. code:: python

   >>> speclist.bin(0.25) # binning our data to a bin width of 0.25 m/z

In FIE-MS null values should concern no more than 3% of the total number
of identified bins. However, imputation is required to streamline the
analysis process (as most multivariate techniques are unable to
accomodate missing data points). To perform value imputation, just use
``value_imputate``:

.. code:: python

   >>> speclist.value_imputate()

Now transforming and normalisating the the spectrum objects in an
samples independent fashion can be done using the following:

.. code:: python

   >>> speclist.transform()
   >>> speclist.normalise()

Once completed, you are now free to export the data to a data matrix:

.. code:: python

   >>> speclist.to_csv("/path/to/proc_metabo.csv", output_type="matrix")

This should give you something akin to:

+-----------+------+-----+------+------+---+
| Sample ID | M0   | M1  | M2   | M3   | … |
+===========+======+=====+======+======+===+
| Sample 1  | 213  | 634 | 3213 | 546  | … |
+-----------+------+-----+------+------+---+
| Sample 2  | 132  | 34  | 713  | 6546 | … |
+-----------+------+-----+------+------+---+
| Sample 3  | 1337 | 42  | 69   | 420  | … |
+-----------+------+-----+------+------+---+
