from dimspy import SpectrumList
from dimspy.Analysis import DataAnalysisObject
from dimspy.Results import ResultsList

import dimspy.Graphics as grfx

import pandas as pd

def get_metadata_dict(fp):
    return pd.read_excel(fp, index_col=0)["FEV1 %  Pred"].dropna()


def load_spectrum_list(fp):
    sl = SpectrumList()
    sl.from_pickle(fp)
    return sl

if __name__ == "__main__":

    mtd = get_metadata_dict("/home/keo7/Desktop/experiments/denisa_saliva/data/Metadata.xlsx")

    #mtd = mtd[mtd.isin(["HC", "COPD", "LC"])]

    sl = load_spectrum_list("/home/keo7/Desktop/experiments/denisa_saliva/data/DIMSpy/positive/processed.pkl")

    da = DataAnalysisObject(sl)
    rl = ResultsList()

    da.linear_regression(mtd)

    exit(0)

    da.principle_components_analysis(mtd, show=True)



    rl.append(da.anova(mtd))

    limit_values = {
        "anova": {
            "p-value": "< 0.05"
        }
    }

    rl.variable_limiter(da, limit_values)

    exit(0)

    da.principle_components_analysis(mtd, show=True)

    exit(0)


    rl.append(da.t_test(mtd))

    limit_values = {
        "t-test": {
            "p-value": "< 0.05"
        }
    }

    rl.variable_limiter(da, limit_values)


    grfx.box_plots(da, mtd, "./test.pdf")

    exit(0)
    rl.append(da.lda(mtd, cv="loo"))

    limit_values = {
        "LDA (Variables)" : {
            "AUC" : "> 0.7"
        }
    }

    rl.variable_limiter(da, limit_values)


    r = da.lda(mtd, cv="loo", type="all")

    grfx.roc(r, True)

    da.principle_components_analysis(mtd, show=True)