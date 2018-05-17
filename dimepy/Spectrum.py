# -*- coding: utf-8 -*-
import os
from collections import defaultdict
from copy import copy
import numpy as np
import operator
import pymzml
from scipy.sparse import csc_matrix, eye, diags
from scipy.sparse.linalg import spsolve
import warnings
import matplotlib.pyplot as plt
from operator import itemgetter

class Spectrum(object):
    """Summary of class here.

    Longer class information...

    Attributes:
        something: description of said thing.

    """
    polarity_dict = {"POSITIVE": "MS:1000130", "NEGATIVE": "MS:1000129"}

    _loaded = False

    _normalised = False
    _transformed = False
    _scaled = False
    _baseline_corrected = False

    _injection_order = None

    masses = np.array([])
    intensities = np.array([])

    __raw_spectrum = None

    def __init__(self,
                 file_path=None,
                 id=None,
                 polarity=None,
                 type="peaks",
                 min_mz=50.0,
                 max_mz=5000,
                 apex=False,
                 snr_estimator="median",
                 max_snr=2.5,
                 injection_order=None):
        """Inits a Spectrum.

        Forms a Spectrum object, contained masses and intensities..

        Args:
            file_path:
            id:
            polarity:
            type:
            apex:
            injection_order:

        """
        self.file_path = file_path
        if id is None:
            self._get_id_from_fp()
        else:
            self.id = id

        if injection_order is not None:
            self._injection_order = int(injection_order)
        else:
            self._injection_order = 0
        self.apex = apex
        self.max_snr = max_snr
        self.snr_estimator = snr_estimator
        self.type = type
        self.min_mz = min_mz
        self.max_mz = max_mz
        self.polarity = polarity
        if self.file_path != None:
            self._load_from_file()
        self._normalised = False
        self._transformed = False
        self._baseline_corrected = False
        self._scaled = False

    def _get_id_from_fp(self):
        """Description of method.

        Longer description of method.
        """
        self.id = os.path.splitext(os.path.basename(self.file_path))[0]

    # Over-egging the WhittakerSmoothing, need to take a look.


    def align(self, ppm=2.0, edge_extend=10):
        from scipy.cluster import hierarchy
        from scipy.spatial import distance
        # ppm: the hierarchical clustering cutting height, i.e., ppm range for each aligned mz value. Default = 2.0
        # ppm range for the edge blocks. Default = 10

        # multiprocess cluster
        def _cluster_peaks_mp(params):
            return _cluster_peaks(*params)


        # single cluster
        def _cluster_peaks(mzs, ppm, distype='euclidean', linkmode='centroid'):
            if len(mzs) == 0:
                return np.array([])
            if len(mzs) == 1:
                return np.zeros_like(mzs, dtype=int).reshape((-1, 1))

            m = np.column_stack([mzs])
            mdist = distance.pdist(m, metric=distype)
            outer_mzs = np.add.outer(mzs, mzs)
            np.fill_diagonal(outer_mzs, 0)
            avg_mz_pair = np.divide(outer_mzs, 2)
            mdist_mz_pair = distance.squareform(avg_mz_pair)
            relative_errors = np.multiply(mdist_mz_pair, 1e-6)

            with np.errstate(divide='ignore', invalid='ignore'):  # using errstate context to avoid seterr side effects
                m_mass_tol = np.divide(mdist, relative_errors)
                m_mass_tol[np.isnan(m_mass_tol)] = 0.0
            z = hierarchy.linkage(m_mass_tol)

            # cut tree at ppm threshold & order matches the order of mzs
            return hierarchy.cut_tree(z, height=ppm)

        def _cluster_peaks_map(block_size=1000, fixed_block=True, ncpus=1):
            mzs = np.array(self.masses)

            if not np.all(mzs[1:] >= mzs[:-1]):
                raise ValueError('mz values not in ascending order')
            if not 1 <= block_size <= len(mzs):
                # logging.warning('block size (%d) not in range [1, #peaks (%d)]' % (block_size, len(mzs)))
                block_size = min(max(block_size, 1), len(mzs))

            # split blocks
            if fixed_block:
                sids = range(block_size, len(mzs), block_size)
            else:
                umzs, urids = np.unique(mzs, return_index=True)
                sids = [urids[i] for i in range(block_size, len(umzs), block_size)]  # appx size
            if len(sids) > 0 and sids[-1] == len(mzs) - 1:
                sids = sids[:-1]  # no need to cluster the final peak separately

            # create mapping pool
            def _mmap(f, p):
                ppool = Pool(cpu_count() - 1 if ncpus is None else ncpus)
                rets = ppool.map(f, p)
                ppool.close()  # close after parallel finished
                return rets

            def _smap(f, p):
                return map(f, p)

            def _pmap(f, p):
                largechk = filter(lambda x: len(x[0]) > 1E+5, p)
                if len(largechk) > 0:
                    raise RuntimeError('Some of the clustering chunks contain too many peaks: \n%s' %
                        join(map(lambda x: 'mz range [%.5f - %.5f] ... [%d] peaks' % (min(x[0]),max(x[0]),len(x[0])), largechk), '\n'))
                return (_smap if ncpus == 1 or cpu_count() <= 2 else _mmap)(f, p)

            # align edges
            eeppm = edge_extend * ppm * 1e-6

            _rng = lambda x: (lambda v: np.where(np.logical_and(x - v < mzs, mzs <= x + v))[0])(eeppm * x)
            erngs = [_rng(mzs[i]) for i in sids]
            overlap = [p[-1] >= s[0] for p, s in zip(erngs[:-1], erngs[1:])]  # in case edges have overlap
            if True in overlap:
                logging.warning('[%d] edge blocks overlapped, consider increasing the block size' % (sum(overlap) + 1))
                erngs = reduce(lambda x, y: (x[:-1] + [np.unique(np.hstack((x[-1], y[0])))]) if y[1] else x + [y[0]],
                               zip(erngs[1:], overlap), [erngs[0]])
                sids = [sids[0]] + [s for s, o in zip(sids[1:], overlap) if not o]

            _cids = _pmap(_cluster_peaks_mp, [(mzs[r], ppm) for r in erngs])
            eblks = [r[c == c[r == s]] for s, r, c in zip(sids, erngs, map(lambda x: x.flatten(), _cids))]
            ecids = map(lambda x: np.zeros_like(x).reshape((-1, 1)), eblks)

            # align blocks
            # keep () in reduce in case eblks is empty
            brngs = np.array(
                (0,) + reduce(lambda x, y: x + y, map(lambda x: (x[0], x[-1] + 1), eblks), ()) + (len(mzs),)
            ).reshape((-1, 2))

            # in case edges have reached mz bounds
            bkmzs = [mzs[slice(*r)] for r in brngs]
            slimbk = map(lambda x: len(x) == 0 or abs(x[-1] - x[0]) / x[0] < eeppm * 10, bkmzs)
            if np.sum(slimbk) > 0:
                pbrngs = [map(lambda x: min(x, len(mzs) - 1), (r[0], r[-1] - 1 if r[-1] != r[0] else r[-1])) for r in brngs]
                pblns = ['block %d' % i + ': [%f, %f]' % itemgetter(*r)(mzs)
                         for i, (s, r) in enumerate(zip(slimbk, pbrngs)) if s]
                logging.warning('[%d] empty / slim clustering block(s) found, consider increasing the block size\n%s' %
                                (np.sum(slimbk), join(pblns, '\n')))
            bcids = _pmap(_cluster_peaks_mp, [(m, ppm) for m in bkmzs])

            # combine
            cids = [None] * (len(bcids) + len(ecids))
            cids[::2], cids[1::2] = bcids, ecids
            return cids

        def _cluster_peaks_reduce(clusters):
            return reduce(lambda x, y: np.vstack((x, y + np.max(x) + 1)), filter(lambda x: len(x) > 0, clusters)).flatten()


        def _align_peaks(cids, pids):
            pids = np.array(pids)
            # encode string id list to continuous values and search uniques
            def _idsmap(ids):
                _, ri, vi = np.unique(ids, return_index=True, return_inverse=True)
                sri = np.argsort(ri)  # ensure order
                return np.argsort(sri)[vi], ids[ri[sri]]

            (mcids, mpids), (ucids, upids) = zip(*map(_idsmap, (cids, pids)))

            print len(upids)
            print len(ucids)
            # count how many peaks from same sample being clustered into one peak
            cM = np.zeros([99267, 54202])
            for pos, count in Counter(zip(mpids, mcids)).items():
                cM[pos] = count

            # fill all the attributes into matrix
            def _avg_am(a):
                aM = np.zeros(map(len, (upids, ucids)))
                for p, v in zip(zip(mpids, mcids), a):
                    aM[p] += v
                with np.errstate(divide='ignore', invalid='ignore'):
                    aM /= cM
                aM[np.isnan(aM)] = 0
                return aM

            def _cat_am(a):
                aM = [[[] for _ in ucids] for _ in upids]
                for (r, c), v in zip(zip(mpids, mcids), a):
                    aM[r][c] += [str(v)]
                aM = [[join(val, ',') for val in ln] for ln in aM]
                return np.array(aM)

            def _fillam(a):
                alg = _avg_am if a.dtype.kind in ('i', 'u', 'f') else \
                      _cat_am if a.dtype.kind in ('?', 'b', 'a', 'S', 'U') else \
                      lambda x: logging.warning('undefined alignment behaviour for [%s] dtype data') # returns None
                return alg(a)

            attrMs = map(_fillam, attrs)

            # sort mz values, ensure mzs matrix be the first
            sortids = np.argsort(np.average(attrMs[0], axis=0, weights=attrMs[0].astype(bool)))
            return upids, map(lambda x: x[:, sortids], attrMs + [cM])

        clusters = _cluster_peaks_map()
        cids = _cluster_peaks_reduce(clusters)
        a_pids, a_attrms = _align_peaks(cids, self.masses)

        print a_pids

    def baseline_correction(self, inplace=True, max_iterations=2, lambda_=100):
        """Description of method.

        Longer description of method.
        """
        warnings.warn("This is currently in development...")

        def _WhittakerSmooth(intensities_copy, ones):
            intensities_copy = np.matrix(intensities_copy)
            diag_eye = eye(intensities_copy.size, format="csc")
            diag = diag_eye[1:] - diag_eye[:-1]
            sparse = diags(
                ones, 0, shape=(intensities_copy.size, intensities_copy.size))

            csc_A = csc_matrix(sparse + (lambda_ * diag.T * diag))
            csc_B = csc_matrix(sparse * intensities_copy.T)
            return np.array(spsolve(csc_A, csc_B))

        def _AirPLS():
            intensities_copy = self.intensities
            ones = np.ones(self.intensities.shape[0])
            for index in range(0, max_iterations):
                whittaker_smoothed = _WhittakerSmooth(intensities_copy, ones)
                smoothed_intensities = (intensities_copy - whittaker_smoothed)
                smoothed_sum = np.abs(
                    smoothed_intensities[smoothed_intensities < 0].sum())
                if (smoothed_sum < 0.001 * (abs(intensities_copy)).sum()
                        or index == max_iterations):
                    break
                ones[smoothed_intensities >= 0] = [0]
                ones[smoothed_intensities < 0] = np.exp(
                    index * (smoothed_intensities[smoothed_intensities < 0]) /
                    smoothed_sum)
                ones[0] = np.exp(
                    index *
                    (smoothed_intensities[smoothed_intensities < 0]).max() /
                    smoothed_sum)
                ones[-1] = ones[0]
            return smoothed_intensities

        if self._baseline_corrected is True:
            warnings.warn(
                "It seems like this spectrum has already been baseline corrected!"
            )

        calculated_baseline = _AirPLS()

        baseline_corrected_intensities = []
        baseline_corrected_masses = []

        for index, intensity in enumerate(
                self.intensities - calculated_baseline):
            if intensity > 0:
                baseline_corrected_intensities.append(intensity)
                baseline_corrected_masses.append(self.masses[index])

        baseline_corrected_intensities = np.array(
            baseline_corrected4_intensities)
        baseline_corrected_masses = np.array(baseline_corrected_masses)

        if inplace is True:
            self.intensities = baseline_corrected_intensities
            self.masses = baseline_corrected_masses
            self._baseline_corrected = True
        else:
            return baseline_corrected_masses, baseline_corrected_intensities

    def _normalise(self, method="tic"):
        '''

        :param method:
        :param inplace:
        :return:
        '''

        if self._normalised is False:
            if method.upper() == "TIC":
                sum_intensity = np.nansum(self.intensities)
                normalised_intensities = np.array([(
                    x / sum_intensity) for x in self.intensities]) * 1000
            elif method.upper() == "MEDIAN":
                median_intensity = np.nanmedian(self.intensities)
                normalised_intensities = np.array(
                    [x - median_intensity for x in self.intensities]) * 1000
            else:
                normalised_intensities = self.intensities
        else:
            warnings.warn("Warning: %s already normalised, ignoring" % self.id)
            normalised_intensities = self.intensities

        self.intensities = normalised_intensities
        self._normalised = True

    def _transform(self, method="log10", inplace=True):
        '''

        :param method:
        :param inplace:
        :return:
        '''

        if self._transformed is False:
            if method.upper() == "LOG10":
                transformed_intensities = np.log10(self.intensities)
            elif method.upper() == "CUBE":
                transformed_intensities = np.array(
                    [i**(1. / 3) for i in self.intensities])
            elif method.upper() == "NLOG":
                transformed_intensities = np.log(self.intensities)
            elif method.upper() == "LOG2":
                transformed_intensities = np.log2(self.intensities)
            elif method.upper() == "GLOG":

                def __lognorm(x, min):
                    return np.log2((x + np.sqrt(x**2 + min**2)) / 2)

                transformed_intensities = __lognorm(self.intensities,
                                                    min(self.intensities) / 10)
            else:
                pass
        else:
            warnings.warn("Warning: %s already normalised, ignoring" % self.id)
            transformed_intensities = self.intensities

        if inplace is True:
            self.intensities = transformed_intensities
            self._transformed = True
        else:
            return transformed_intensities

    def _load_from_file(self):
        def __from_mzml():
            def __gen_reader():
                return pymzml.run.Reader(
                    self.file_path,
                    extraAccessions=[["MS:1000129", ["value"]],
                                     ["MS:1000130", ["value"]]])

            def __get_scans_of_interest():
                reader = __gen_reader()
                scans = []
                for scan_number, scan in enumerate(reader):
                    if self.polarity != None:
                        if scan.get(self.polarity_dict[self.polarity.upper()]
                                    ) != None:
                            scans.append(scan_number)
                    else:
                        scans.append(scan_number)

                reader = __gen_reader()
                if self.apex == True:
                    tics = []
                    for scan_number, scan in enumerate(reader):
                        tic = sum(zip(*scan.peaks)[1])
                        tics.append([scan_number, tic])
                    mad = np.mean(np.absolute(zip(*tics)[0] - np.mean(zip(*tics)[0])))
                    scans = [tics[i][0] for i, x in enumerate(zip(*tics)[1]) if x > mad]
                return scans


            def __get_scan(scan_range):
                masses = []
                intensities = []
                reader = __gen_reader()
                for scan_number, scan in enumerate(reader):
                    if scan["ms level"] != None:
                        m, ints = zip(*scan.peaks)
                        # TODO: This but better.
                        nm, nints = [], []
                        for indx, mass in enumerate(m):
                            if mass > self.min_mz and mass < self.max_mz:
                                nm.append(mass)
                                nints.append(ints[indx])
                        m, ints = nm, nints
                        if len(m) > 0:
                            if self.snr_estimator != None:
                                sn_r = np.divide(ints, scan.estimatedNoiseLevel(mode=self.snr_estimator))
                                gq = [i for i, e in enumerate(sn_r) if e < self.max_snr]
                                m = np.array(m)[gq]
                                ints = np.array(ints)[gq]
                            masses.extend(m)
                            intensities.extend(ints)

                masses, intensities = zip(*sorted(zip(masses, intensities)))
                return masses, intensities

            scan_range = __get_scans_of_interest()

            self.masses, self.intensities = __get_scan(scan_range)


        if self.file_path.upper().endswith("MZML"):
            __from_mzml()
        elif self.file_path.upper().endswith("CSV"):
            print "CSV not implemented"
        elif self.file_path.upper().endswith("PKL"):
            print "PKL not implemented"
        else:
            raise Exception()

    def plot(self, show=True, xlim=[], scaled=False, file_path=None):
        plt.figure()
        plt.title(self.id)
        plt.xlabel("Mass-to-ion (m/z)")
        plt.ylim(0, max(self.intensities))
        plt.ylabel("Intensity")

        if xlim == []:
            xlim = [min(self.masses), max(self.masses)]
            plt.ylim(0, max(self.intensities))
        plt.xlim(xlim)
        if scaled is False:
            plt.ylim(0, max(self.intensities))
            plt.plot(self.masses, self.intensities)
            plt.ylabel("Intensity")
        else:
            scaled_intensities = self.normalise(method="tic", inplace=False)
            plt.plot(self.masses, scaled_intensities)
            plt.ylim(0, max(scaled_intensities))
            plt.ylabel("Scaled Intensity")
        plt.tight_layout()
        if file_path is not None:
            plt.savefig(file_path)
        else:
            plt.show()
        plt.clf()
