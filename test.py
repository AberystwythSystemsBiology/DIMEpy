from dimspy import SpectrumList
from dimspy.Analysis import DataAnalysisObject
from dimspy.Results import ResultsList

import pandas as pd

def get_metadata_dict(fp):
    return pd.read_excel(fp, index_col=0)["Category Diagnosis"].dropna()


def load_spectrum_list(fp):
    sl = SpectrumList()
    sl.from_pickle(fp)
    return sl

if __name__ == "__main__":

    mtd = get_metadata_dict("/home/keo7/Desktop/experiments/denisa_saliva/data/Metadata.xlsx")

    mtd = mtd[mtd.isin(["HC", "COPD"])]

    sl = load_spectrum_list("/home/keo7/Desktop/experiments/denisa_saliva/data/DIMSpy/positive/processed.pkl")

    da = DataAnalysisObject(sl)

    rl = ResultsList()

    rl.append(da.t_test(mtd))

    limit_values = {
        "t-test": {
            "p-value": "< 0.05"
        }
    }

    rl.variable_limiter(da, limit_values)

    da.lda(mtd, type="variable", cv="loo")


    #da.pca(mtd, show=True)
