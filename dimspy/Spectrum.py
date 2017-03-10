import numpy as np, operator, pymzml, pickle as pkl, pandas as pd

class Spectrum(object):

    def __init__(self, id, masses=np.array([]), intensities=np.array([])):
        self.id = id
        self.masses = np.array(masses)
        self.intensities = np.array(intensities)

    def pickle(self, fp):
        with open(fp, "wb") as output:
            pkl.dump(self, output, pkl.HIGHEST_PROTOCOL)

    def from_pickle(self, fp):
        with open(fp, "rb") as infile:
            o = pkl.load(infile)
        self.id = o.id
        self.masses = o.masses
        self.intensities = o.intensities

    def to_csv(self, fp="/tmp/spectrum.csv", delimiter=","):
        output = []
        output.append([self.id] + self.masses.tolist())
        output.append([" "] + self.intensities.tolist())
        pd.DataFrame(output).transpose().to_csv(fp, delimiter=delimiter, header=False, index=False)

    def from_csv(self, fp="/tmp/spectrum.csv", delimiter=","):
        data = pd.DataFrame.from_csv(fp).T
        self.masses = data.columns.values
        self.intensities = data.values[0]

    def smooth(self, sigma=3):
        from scipy.ndimage.filters import gaussian_filter1d
        self.intensities = gaussian_filter1d(self.intensities, sigma)

    def correct_baseline(self, lambda_=100, porder=1, max_iterations=15):
        from scipy.sparse import csc_matrix, eye, diags
        from scipy.sparse.linalg import spsolve

        def _WhittakerSmooth(x, w, lambda_, differences=1):
            X = np.matrix(x)
            m = X.size
            i = np.arange(0, m)
            E = eye(m, format='csc')
            D = E[1:] - E[:-1]
            W = diags(w, 0, shape=(m, m))
            A = csc_matrix(W + (lambda_ * D.T * D))
            B = csc_matrix(W * X.T)
            background = spsolve(A, B)
            return np.array(background)

        def _AirPLS():
            m = self.intensities.shape[0]
            w = np.ones(m)
            for i in range(1, max_iterations + 1):
                z = _WhittakerSmooth(self.intensities, w, lambda_, porder)
                d = self.intensities - z
                dssn = np.abs(d[d < 0].sum())
                if (dssn < 0.001 * (abs(self.intensities)).sum() or i == max_iterations):
                    break
                w[d >= 0] = 0
                w[d < 0] = np.exp(i * np.abs(d[d < 0]) / dssn)
                w[0] = np.exp(i * (d[d < 0]).max() / dssn)
                w[-1] = w[0]
            return z

        baseline = _AirPLS()

        bc_i = []
        bc_m = []
        for index, intensity in enumerate(self.intensities - baseline):
            if intensity > 0:
                bc_i.append(self.intensities[index])
                bc_m.append(self.masses[index])

        self.intensities = np.array(bc_i)
        self.masses = np.array(bc_m)


    def normalise(self, method="tic"):
        if method == "tic":
            sum_intensity = np.nansum(self.intensities)
            median_intensity = np.nanmedian(self.intensities)
            self.intensities = np.array([(x / sum_intensity) * median_intensity for x in self.intensities])
        else:
            pass

    def transform(self, method="log10"):
        if method == "log10":
            self.intensities = np.log10(self.intensities)

    def get_peaks(self, delta=3):
        def _detect_peaks():
            # TODO: Legacy code - badly needs a rewrite.
            maxtab = []
            mintab = []
            x = range(len(self.intensities))
            v = np.asarray(self.intensities)
            mn, mx = np.Inf, -np.Inf
            mnpos, mxpos = np.NaN, np.NaN
            lookformax = True
            for i in np.arange(len(v)):
                this = v[i]
                if this > mx:
                    mx = this
                    mxpos = x[i]
                if this < mn:
                    mn = this
                    mnpos = x[i]
                if lookformax:
                    if this < mx - delta:
                        maxtab.append(mxpos)
                        mn = this
                        mnpos = x[i]
                        lookformax = False
                else:
                    if this > mn + delta:
                        mintab.append(mnpos)
                        mx = this
                        mxpos = x[i]
                        lookformax = True

            return np.array(maxtab)

        peaks = _detect_peaks()
        self.masses = self.masses[peaks]
        self.intensities = self.intensities[peaks]

    def from_mzml(self, filepath, polarity=None, scan_range="all", peak_type="peaks", ms1_p=5e-6, msn_p=5e-6):

        def _get_polarity():
            reader = pymzml.run.Reader(filepath, MS1_Precision=ms1_p, MSn_Precision=msn_p,
                                       extraAccessions=[('MS:1000129', ['value']), ('MS:1000130', ['value'])])

            polarity_dict = {"negative" : "MS:1000130", "positive" : "MS:1000129"}
            scans_of_interest = []
            for scan_number, scan in enumerate(reader):
                if scan.get(polarity_dict[polarity]) != None:
                    scans_of_interest.append(scan_number)
            return scans_of_interest

        def _return_apex():
            reader = pymzml.run.Reader(filepath, MS1_Precision=ms1_p, MSn_Precision=msn_p)
            tic_scans = []
            for scan_number, scan in enumerate(reader):
                if scan_number in scans_of_interest:
                    tic_scans.append([scan["total ion current"], scan_number])
            tics = [x[0] for x in tic_scans]
            mad = np.mean(np.absolute(tics - np.mean(tics))) * 3
            scan_range = [x[1] for x in tic_scans if x[0] > mad]
            return scan_range

        def _create_spectrum():
            reader = pymzml.run.Reader(filepath, MS1_Precision=ms1_p, MSn_Precision=msn_p)
            sample_spectrum = pymzml.spec.Spectrum(measuredPrecision=msn_p)
            for scan_number, scan in enumerate(reader):
                if scan_number in scans_of_interest:
                    sample_spectrum += scan
            return sample_spectrum

        def _pymzl_spectrum_to_spectrum():
            if peak_type == "centroided":
                spectrum = [[masses, intensities] for masses, intensities in sample_spectrum.centroidedPeaks]
            elif peak_type == "reprofiled":
                spectrum = [[masses, intensities] for masses, intensities in sample_spectrum.reprofiledPeaks]
            else:
                spectrum = [[masses, intensities] for masses, intensities in sample_spectrum.peaks]
            spectrum = sorted(spectrum, key=operator.itemgetter(0))

            masses = np.array([x[0] for x in spectrum])
            intensities = np.array([x[1] for x in spectrum])

            self.masses = masses
            self.intensities = intensities

        scans_of_interest = _get_polarity()
        if scan_range == "apex":
            scans_of_interest = _return_apex()

        sample_spectrum = _create_spectrum()
        _pymzl_spectrum_to_spectrum()