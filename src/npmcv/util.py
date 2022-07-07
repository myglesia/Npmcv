from itertools import zip_longest
from scipy import stats
import numpy as np


# Utility and stats functions


def grouper(iterable, n, fillvalue=None):
    '''Collect data into fixed-length chunks or blocks
    grouper('ABCDEFG', 3, 'x') --> ABC DEF Gxx'''
    # Adapted from: https://docs.python.org/3/library/itertools.html
    args = [iter(iterable)] * n
    return zip_longest(*args, fillvalue=fillvalue)


def remove_outliers(results):
    '''results:    {'img': [CVs]}'''
    clean_results = {}
    outliers = {}
    for key, val in results.items():
        val = [x for x in val if not np.isnan(x)]  # remove 'NaN'
        if val:
            clean, out = grubbs(val)
            clean_results[key] = clean
            outliers[key] = out

    return clean_results, outliers


def grubbs(X, alpha=0.05):
    '''
    Performs Grubbs' test for outliers recursively until the null hypothesis is
    true.
    Parameters
    ----------
    X : ndarray
        A numpy array to be tested for outliers.
    alpha : float
        The significance level.
    Returns
    -------
    X : ndarray
        The original array with outliers removed.
    outliers : ndarray
        An array of outliers.
    '''
    Z = stats.zscore(X, ddof=1)  # returns ndarray of Z-score
    N = len(X)

    def G(Z): return np.abs(Z).argmax()  # returns indices of max values
    def t_crit(N): return stats.t.isf(alpha / (2. * N), N - 2)

    def G_crit(N): return (N - 1.) / np.sqrt(N) * \
        np.sqrt(t_crit(N)**2 / (N - 2 + t_crit(N)**2))

    outliers = np.array([])
    while np.amax(np.abs(Z)) > G_crit(N):
        outliers = np.r_[outliers, X[G(Z)]]
        # remove outlier from array
        X = np.delete(X, G(Z))
        # repeat Z score
        Z = stats.zscore(X, ddof=1)
        N = len(X)

    return X, outliers