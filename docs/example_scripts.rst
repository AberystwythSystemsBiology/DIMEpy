Example Scripts
===============

Here's a quick run through of key functionality in DIMEpy.l

.. code:: python

    import dimepy
    import os

    data_dir = "/path/to/mzMLs"

    sl = dimepy.SpectrumList()

    for fn in os.listdir(data_dir):
         # Appending the file name to the data directory
         fp = os.path.join(data_dir,fn)

         # We only have two classes, denoted by the first letter of the file
         strat = fn[0]

         # Note how I'm applying mean-based snr estimation
         spec = dimepy.Spectrum(fp, fn, stratification=strat, snr_estimator="mean")

         # As these are polarity-switching, I'm limiting to positive
         spec.limit_polarity("positive")

         # A threshold of 3 seemed fine from earlier.
         spec.limit_infusion(3)

         # Make sure I load the scans
         spec.load_scans()

         # Appending the Spectrum object to the SL
         sl.append(spec)

    sl.detect_outliers(3, verbose=True)

    sl.bin(0.25)
    sl.value_imputate()
    sl.transform()
    sl.normalise()

    sl.to_csv("/path/to/output.csv", output_type="matrix")