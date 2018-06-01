'''
    def value_imputation(self, method="basic", threshold=0.5, inplace=True):
        def _remove_bins_by_threshold():
            df = self._spectrum_list.flatten_to_dataframe()

            nc = df.isnull().sum()
            _threshold = len(df.index.values) * threshold
            nc = nc[nc <= _threshold]
            return df[nc.index.values]

        def _value_imputation(df):
            if method.upper() == "KNN":
                imp = Imputer(axis=1)
                imputated = imp.fit_transform(df)
                df = pd.DataFrame(
                    imputated, columns=df.columns, index=df.index)
            else:
                for sample in df.index:
                    intensities = df.ix[sample].values
                    if method.upper() == "BASIC":
                        filler = np.nanmin(intensities) / 2
                    elif method.upper() == "MEAN":
                        filler = np.nanmean(intensities)
                    elif method.upper() == "MIN":
                        filler = np.nanmin(intensities)
                    elif method.upper() == "MEDIAN":
                        filler = np.nanmedian(intensities)
                    df.ix[sample] = df.ix[sample].replace(np.nan, filler)

            return df

        if method.upper() == "ALL":
            threshold = 0
            df = _remove_bins_by_threshold()
        else:
            df = _remove_bins_by_threshold()
            df = _value_imputation(df)

        if inplace is True:
            self.pandas_to_spectrum(df)
            self._value_imputated = True
        else:
            return df
'''
