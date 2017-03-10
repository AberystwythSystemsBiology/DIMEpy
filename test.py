from dimspy import Spectrum, SpectrumList, SpectrumListProcessor

file_path = "/home/keo7/Desktop/experiments/denisa_saliva/data/mzMLs/"
scan_range = "apex"
type = "peaks"

spec_list = SpectrumList()

spec = Spectrum(id="A002")
spec.from_mzml(file_path+"0001-A002-160824-a.mzML", "negative", scan_range, type)
spec_list.append(spec)

spec2 = Spectrum(id="Ctrl01")
spec2.from_mzml(file_path+"Ctrl01.mzML", "negative", scan_range, type)
spec_list.append(spec2)

sl_processor = SpectrumListProcessor(spec_list)
sl_processor.binning(0.125)
sl_processor.value_imputation(threshold=0.5)

processed_spec_list = sl_processor.to_spectrumlist()

processed_spec_list.to_excel()