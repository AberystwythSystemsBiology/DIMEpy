from dimspy import Spectrum, SpectrumList, SpectrumListProcessor

import matplotlib.pyplot as plt

file_path = "/home/keo7/Desktop/experiments/denisa_saliva/data/mzMLs/"
scan_range = "apex"
type = "peaks"

spec_list = SpectrumList()

spec = Spectrum(id="A002")
spec.from_mzml(file_path+"0001-A002-160824-a.mzML", "negative", scan_range, type)

spec_list.append(spec)

sl_processor = SpectrumListProcessor(spec_list)
sl_processor.outlier_detection()
sl_processor.smooth()
sl_processor.correct_baseline()
sl_processor.binning(0.125)
sl_processor.value_imputation(threshold=0.5)
sl_processor.normalise()
sl_processor.transform()

processed_spec_list = sl_processor.to_spectrumlist()

plt.figure()
for spectrum in processed_spec_list.to_list():
    plt.plot(spectrum.masses, spectrum.intensities, label=spectrum.id)
plt.legend()
plt.show()