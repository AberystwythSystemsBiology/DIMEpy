from dimspy import Spectrum, SpectrumList, SpectrumListProcessor, Analysis

spectrum_list = SpectrumList()

spectrum = Spectrum(id="spectrum_one",
                    masses=[10.0, 10.1, 10.2, 10.3, 10.4, 10.5, 10.6, 10.7, 10.8, 10.9, 11.0],
                    intensities=[30, 10, 15, 10, 15, 20, 25, 10, 5, 10, 100])


spectrum_list.append(spectrum)

spectrum_two = Spectrum(id="spectrum_two",
                        intensities=[9, 10, 13, 14, 20, 100, 23, 10, 50, 90, 10, 20, 30],
                        masses=[9, 9.6, 9.7, 9.8, 9.9, 10.0, 10.1, 10.2, 10.3, 10.4, 10.5, 10.6, 10.7])

spectrum_list.append(spectrum_two)
Analysis.DataAnalysisObject
spectrum_list_processor = SpectrumListProcessor(spectrum_list)

spectrum_list_processor.binning(0.15, statistic="median")
spectrum_list_processor.value_imputation(method="basic", threshold=0.5)

spectrum_list = spectrum_list_processor.to_spectrumlist()

spectrum_list.to_excel("/home/keo7/Desktop/test.xlsx")