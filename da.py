from dimspy import SpectrumListProcessor, SpectrumList, Analysis, Results
import pandas as pd, matplotlib.pyplot as plt
sl = SpectrumList()
sl.from_pickle("/home/keo7/Desktop/experiments/denisa_saliva/data/DIMSpy/positive/processed.pkl")

def get_metadata_dict(fp):
    return pd.read_excel(fp, index_col=0)["Category Diagnosis"].dropna()


da = Analysis.DataAnalysisObject(sl)

rl = Results.ResultsList()

mtd = get_metadata_dict("/home/keo7/Desktop/experiments/denisa_saliva/data/Metadata.xlsx")

mtd = mtd[mtd.isin(["HC", "LC"])]

da.principle_components_analysis(mtd, fp="/home/keo7/Desktop/pca_01_mz_bin.png")

rl.append(da.t_test(mtd))

limit_values = {
    "t-test": {
        "p-value": "< 0.05"
    }
}

rl.variable_limiter(da, limit_values)

da.lda(mtd)

da.principle_components_analysis(mtd, fp="/home/keo7/Desktop/pca_02_mz_bin.png")

plt.show()
