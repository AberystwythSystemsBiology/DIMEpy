from dimspy import Spectrum, SpectrumList, SpectrumListProcessor
bd = "/home/keo7/Desktop/experiments/denisa_saliva/data/mzMLs/"

sl = SpectrumList()

for i in ["Ctrl01.mzML", "Ctrl02.mzML", "Ctrl03.mzML", "Ctrl04.mzML", "Ctrl05.mzML"]:
    spec = Spectrum(id=i)
    spec.from_mzml(bd+i, polarity="negative", scan_range="apex")
    sl.append(spec)

sl.pickle("here.pkl")

sl.from_pickle("here.pkl")

slp = SpectrumListProcessor(sl)

slp.correct_baseline()
slp.center()

slp.value_imputation(threshold=0.5)
slp.normalise()
slp.transform()

spectrum_list = slp.to_spectrumlist()

spectrum_list.to_excel("./here.xlsx")