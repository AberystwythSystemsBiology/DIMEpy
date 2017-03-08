import numpy as np, scipy.stats as stats, operator, pymzml, pickle as pkl
from scipy.ndimage.filters import gaussian_filter1d
from scipy.sparse import csc_matrix, eye, diags
from scipy.sparse.linalg import spsolve


class Spectrum(object):
    def __init__(self, identifier, masses=np.array([]), intensities=np.array([])):
        self.identifier = identifier
        self.masses = np.array(masses)
        self.intensities = np.array(intensities)

    def smooth(self, sigma=1):
       self.intensities  = gaussian_filter1d(self.intensities, sigma)

    def correct_baseline(self, lambda_=100, porder=1, max_iterations=15):
        def WhittakerSmooth(x, w, lambda_, differences=1):
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

        def AirPLS():
            m = self.intensities.shape[0]
            w = np.ones(m)
            for i in range(1, max_iterations + 1):
                z = WhittakerSmooth(self.intensities, w, lambda_, porder)
                d = self.intensities - z
                dssn = np.abs(d[d < 0].sum())
                if (dssn < 0.001 * (abs(self.intensities)).sum() or i == max_iterations):
                    break
                w[d >= 0] = 0
                w[d < 0] = np.exp(i * np.abs(d[d < 0]) / dssn)
                w[0] = np.exp(i * (d[d < 0]).max() / dssn)
                w[-1] = w[0]
            return z

        baseline = AirPLS()

        bc_i = []
        bc_m = []
        for index, intensity in enumerate(self.intensities-baseline):
            if intensity > 0:
                bc_i.append(self.intensities[index])
                bc_m.append(self.masses[index])

        self.intensities = np.array(bc_i)
        self.masses = np.array(bc_m)

    def bin(self, bs=0.25, removena=True):
        bins = np.arange(round(min(self.masses)), round(max(self.masses)), step=bs)
        b_intensities, b_masses, b_num = stats.binned_statistic(self.masses, self.intensities, bins=bins)

        binned_masses = []
        binned_intensities = []

        if removena == True:
            for index, intensity in enumerate(b_intensities):
                if np.isnan(intensity) != True:
                    binned_masses.append(b_masses[index])
                    binned_intensities.append(b_intensities[index])
        self.masses = np.array(binned_masses)
        self.intensities = np.array(binned_intensities)

    def normalise(self, method="tic"):
        if method == "tic":
            sum_intensity = np.sum(self.intensities)
            median_intensity = np.median(self.intensities)
            self.intensities = np.array([(x / sum_intensity) * median_intensity for x in self.intensities])
        else:
            pass

    def transform(self, method="log10"):
        if method == "log10":
            self.intensities = np.log10(self.intensities)

    def value_imputation(self):
        intensities = self.intensities
        intensities[np.isnan(intensities)] = (np.nanmin(intensities)/2)
        self.intensities = intensities


    def pickle(self, fp):
        with open(fp, "wb") as output:
            pkl.dump(self, output, pkl.HIGHEST_PROTOCOL)

    def from_pickle(self, fp):
        with open(fp, "rb") as infile:
            o = pkl.load(infile)
        self.identifier = o.identifier
        self.masses = o.masses
        self.intensities = o.intensities

    def from_mzml(self, filepath, polarity=None, scan_range="all", peak_type="peaks", ms1_p=5e-6, msn_p=5e-6):

        def get_apex():
            tic_scans = []
            for scan_number, spectrum in enumerate(reader):
                try:
                    if p_d[polarity] in spectrum["filter string"]:
                        tic_scans.append([spectrum["total ion current"], scan_number])
                except KeyError:
                    pass
            tics = [x[0] for x in tic_scans]
            mad = np.mean(np.absolute(tics - np.mean(tics))) * 3
            scan_range = [x[1] for x in tic_scans if x[0] > mad]
            return scan_range

        p_d = {"positive": "+ p",
               "negative": "- p"}
        if polarity not in p_d.keys():
            return
        reader = pymzml.run.Reader(filepath, MS1_Precision=ms1_p, MSn_Precision=msn_p)
        if scan_range == "all":
            scan_range = []
            for scan_number, spectrum in enumerate(reader):
                try:
                    if p_d[polarity] in spectrum["filter string"]:
                        scan_range.append(scan_number)
                except KeyError:
                    pass
        elif scan_range == "apex":
            scan_range = get_apex()
        elif type(scan_range) == list:
            pass

        reader = pymzml.run.Reader(filepath, MS1_Precision=ms1_p, MSn_Precision=msn_p)

        sample_spectrum = pymzml.spec.Spectrum(measuredPrecision=msn_p)

        for scan_number, scan_spectrum in enumerate(reader):
            if scan_number in scan_range:
                sample_spectrum += scan_spectrum

        if peak_type == "peaks":
            spectrum = [[masses, intensities] for masses, intensities in sample_spectrum.peaks]
        elif peak_type == "centroided":
            spectrum = [[masses, intensities] for masses, intensities in sample_spectrum.centroidedPeaks]
        elif peak_type == "reprofiled":
            spectrum = [[masses, intensities] for masses, intensities in sample_spectrum.reprofiledPeaks]

        spectrum = sorted(spectrum, key=operator.itemgetter(0))

        masses = np.array([x[0] for x in spectrum])
        intensities = np.array([x[1] for x in spectrum])

        self.masses = masses
        self.intensities = intensities