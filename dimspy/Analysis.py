import pandas as pd, numpy as np, dimspy_utils
from Results import Result
from sklearn.metrics import roc_auc_score, accuracy_score

class DataAnalysisObject(object):
    def __init__(self, SpectrumList):
        data_frame = []
        for spectrum in SpectrumList.spectrum_list:
            df = pd.DataFrame(spectrum.intensities).T
            df.columns = spectrum.masses
            df.index  = [spectrum.id]
            data_frame.append(df)
        data_frame = pd.concat(data_frame, axis=0)
        self.data_frame = data_frame


    def _append_class(self, class_df):
        return pd.concat([class_df, self.data_frame], axis=1).dropna()

    def principle_components_analysis(self, class_df=None, fp=None, show=False):
        from sklearn.decomposition import PCA
        import matplotlib.pyplot as plt

        def _apply_pca():
            pca = PCA(whiten=True)
            X_r = pca.fit_transform(data)
            return X_r, pca

        def _plot_pca():
            fig, axis = plt.subplots()
            for i, (c, target_name) in enumerate(zip("rgbykm", set(labels))):
                plt.scatter(X_r[labels == target_name, 0],
                            X_r[labels == target_name, 1],
                            color=c, label=target_name, alpha=0.9, s=25)
                plt.xlabel("PC 1 (% .1f %%)" % (
                    pca.explained_variance_ratio_[0] * 100.0))
                plt.ylabel("PC 2 (% .1f %%)" % (
                    pca.explained_variance_ratio_[1] * 100.0))
                for n, x, y in zip(
                        (labels == target_name).nonzero()[0],
                        X_r[labels == target_name, 0],
                        X_r[labels == target_name, 1]):
                    axis.annotate(ids[n], (x, y), va="top", ha="left", fontsize='small')
            if class_df is not None:
                plt.legend(numpoints=1)
            if fp != None:
                plt.savefig(fp)
            if show != False:
                plt.show()

        if class_df is None:
            df = self.data_frame
            data = df[df.columns[0:]].values
            labels = np.ones(len(df.index))
            ids = df.index.values
        else:
            df = self._append_class(class_df)
            data = df[df.columns[1:]].values
            labels = df[df.columns[0]].values
            ids = df.index.values

        X_r, pca = _apply_pca()
        _plot_pca()

    def t_test(self, class_df):
        from scipy.stats import ttest_ind
        classes = list(set(class_df.values))

        if len(classes) > 2:
            return

        df = self._append_class(class_df)
        a_all = df[df[df.columns[0]].isin([classes[0]])]
        b_all = df[df[df.columns[0]].isin([classes[1]])]

        results = []

        for variable in df.columns[1:]:
            a = a_all[variable].values
            b = b_all[variable].values
            t_stat, p_value = ttest_ind(a, b, equal_var=False)
            results.append([variable, p_value, t_stat])

        results = pd.DataFrame([x[1:] for x in results],
                           index=[x[0] for x in results],
                           columns=["p-value", "t-statistic"])

        return Result(results, "t-test")


    def anova(self, class_df):
        pass


    def lda(self, class_df, type="variable", cv=None):
        from sklearn.discriminant_analysis import LinearDiscriminantAnalysis

        df = self._append_class(class_df)

        def _run_lda(data, labels, train=None, test=None):
            clf = LinearDiscriminantAnalysis()
            if train is None:
                clf.fit(data, labels)
                prediction = clf.predict(data)
                score = clf.decision_function(data)
            else:
                clf.fit(data[train], labels[train])
                prediction = clf.predict(data[test])
                score = clf.decision_function(data[test])
            return prediction, score

        def _variables():
            results = []
            data, labels = dimspy_utils._prep_data(df, labels=True, binarise=True)

            for variable in df.columns[1:]:
                y_true = []
                y_predict = []
                y_score = []
                data = df[variable].values.reshape(len(labels), 1)
                if cv != None:
                    c_validator = dimspy_utils._cross_validation(data, cv)
                    for train, test in c_validator:
                        p, s = _run_lda(data, labels, train, test)
                        y_true.extend(labels[test])
                        y_predict.extend(p)
                        y_score.extend(s)
                else:
                    p, s = _run_lda(data, labels)
                    y_true.extend(labels)
                    y_predict.extend(p)
                    y_score.extend(s)

                auc = roc_auc_score(y_true, y_score)
                accuracy = accuracy_score(y_true, y_predict)
                results.append([variable, auc, accuracy])
            results = pd.DataFrame([x[1:] for x in results], index=[x[0] for x in results], columns=["AUC", "Accuracy"])
            return Result(results, "LDA (Variables)")

        def _all():
            y_true = []
            y_predict = []
            y_score = []
            data, labels = dimspy_utils._prep_data(df, labels=True, binarise=True)
            if cv != None:
                c_validator = dimspy_utils._cross_validation(data, cv)
                for train, test in c_validator:
                    p, s = _run_lda(data, labels, train, test)
                    y_true.extend(labels[test])
                    y_predict.extend(p)
                    y_score.extend(s)
            else:
                p, s = _run_lda(data, labels)
                y_true =[x[0] for x in labels.tolist()]
                y_predict = p
                y_score = s
            auc = roc_auc_score(y_true, y_score)
            accuracy = accuracy_score(y_true, y_predict)
            result = pd.DataFrame([auc, accuracy]).T
            result.columns = ["AUC", "Accuracy"]
            result = Result(result, "LDA (All)")
            result.y_true = y_true
            result.y_pred = y_predict
            result.y_scores = y_score
            return result


        if type == "variable":
            result = _variables()
        elif type == "all":
            result = _all()
        return result