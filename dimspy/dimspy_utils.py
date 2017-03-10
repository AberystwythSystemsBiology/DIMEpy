def _cross_validation(data, type, n_folds=None):
    from sklearn.model_selection import LeaveOneOut, StratifiedKFold
    if type == "loo":
        loo = LeaveOneOut()
        return loo.split(data)
    if type == "skf":
        skf = StratifiedKFold()
        return skf.split(data, n_folds)

def _prep_data(df, labels=False, binarise=False):
    if labels == False:
        return df.values
    else:
        data = df[df.columns[1:]].values
        labels = df[df.columns[0]].values
        if binarise == True:
            from sklearn.preprocessing import LabelBinarizer
            lb = LabelBinarizer()
            labels = lb.fit_transform(labels)
        return data, labels

