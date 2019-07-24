#!/usr/bin/env python3

import sys
import os
import time
import re
from itertools import zip_longest
from collections import defaultdict

import mahotas as mh
import pandas as pd
import numpy as np

from PIL import Image
from scipy import stats, misc
from scipy import ndimage as ndi
from skimage.measure import regionprops
from skimage.filters import gaussian, threshold_li


def main(argv):
    '''npmcv <directory>
    DAPI and NPM1 channels from a microscope image must be split into separate tif files. Files should be listed sequentially.

    e.g. "Image01_dapi.tif", "Image01_npm.tif", "Image02_dapi.tif", "Image02_npm.tif",...
    '''
    t = time.process_time()
    directory = argv

    img_files = sorted([f for f in os.listdir(directory)
                        if f.endswith('.tif')], key=lambda f: f.lower())

    pattern = re.compile(r"(?P<slide>.*)\.lif.*")
    raw_results = defaultdict(list)
    for dapi, npm1 in grouper(img_files, 2):
        img_voc = sip(os.path.join(directory, dapi),
                      os.path.join(directory, npm1))
        r = pattern.match(npm1)
        raw_results[r.group(1)].extend(img_voc)

    results, outliers = remove_outliers(raw_results)


    # handling save files
    file_name = os.path.basename(os.path.abspath(directory))

    def dtdf(d): return pd.DataFrame(dict([(k, pd.Series(v)) for k, v in d.items()]))

    df = dtdf(results)
    df.to_csv('{}.csv'.format(file_name), index=False)
    dtdf(raw_results).to_csv('{}_RAW.csv'.format(file_name), index=False)
    dtdf(outliers).to_csv('{}_OUT.csv'.format(file_name), index=False)

    print('\x1b[1;36m' +
          ' {n} Cells Counted.'.format(n=str(sum(df.count())) + '\x1b[0m'))

    elapsed_time = (time.process_time() - t) / 60
    print('\x1b[1;32m' +
          'Total Runtime: {:.4} mins.'.format(elapsed_time) + '\x1b[0m')


def sip(dapi, npm1):
    '''Single image Processing

    Returns
    -------
        results:    list
            list of stats for one image
    '''
    print('\x1b[4m' + 'Processing: {}'.format(os.path.basename(npm1)) + '\x1b[0m')
    with Image.open(dapi) as d, Image.open(npm1) as n:
        raw_dapi = np.array(Image.open(dapi))
        raw_npm1 = np.array(Image.open(npm1))

        label = cell_segmentation(raw_dapi)

        results = []
        if label is None:
            results.append(np.nan)
        else:
            corr_npm1 = bg_correction(raw_npm1)
            cells = img2cells(label, corr_npm1)
            for i, (cell, mask) in enumerate(cells):
                cv = stats.variation(cell[mask], axis=None)
                results.append(cv)
                if cv > 0.75:
                    #verb(i, cell, mask, npm1)
                    continue

    return results


def verb(i, cell, mask, ppath):
    img = np.zeros(cell.shape, dtype=cell.dtype)
    img[:, :] = mask * cell

    name = os.path.dirname(
        ppath) + '/in_cells/{0}_c{1}.tif'.format(os.path.basename(ppath), i)
    os.makedirs(os.path.dirname(name), exist_ok=True)
    misc.toimage(img).save(name)


def cell_segmentation(img):
    '''Prepares a mask labeling each cell in the DAPI images,
    removing artifacts and cut off cells.

    Returns
    -------
    cleaned:    label
        labeled binary image for a single slide
    '''


    # Li Threshold...
    im = gaussian(img.astype(float), sigma=5)  # time!
    bin_image = (im > threshold_li(im))
    labeled, nr_objects = mh.label(bin_image)

    # NOTE Change these values first
    cleaned, n_left = mh.labeled.filter_labeled(
        labeled, remove_bordering=True, min_size=10000, max_size=25000)

    print("Segmentation Calculation...")

    print('Final Number of cells: ' +
          '\x1b[3;31m' + '{}'.format(n_left) + '\x1b[0m')
    # print('Done!\n')

    if n_left == 0:
        return None
    else:
        return cleaned


def bg_correction(img):
    '''Background correction for NPM1 image
    Subtracts the mean background intensity from the entire image

    Returns
    -------
    corrected:  ndarray

    '''
    im = gaussian(img.astype(float), sigma=10)  # time!
    bin_image = (im > threshold_li(im))
    background_mean = img[bin_image == 0].mean()

    # Type change required for saturation
    corrected = (img.astype(np.float) - int(background_mean)).clip(min=0)

    return corrected


def img2cells(labeled, npm):
    '''Separates cells in NPM1 image based on labeled mask from DAPI image
    Returns
    -------
        (cropped_images, cropped_binary):    zip
    '''
    cropped_img = []
    cropped_bin = []

    for region in regionprops(labeled):
        minr, minc, maxr, maxc = region.bbox
        m = labeled[:, :] == region.label
        mask = ndi.binary_fill_holes(m)
        cropped_img.append(npm[minr:maxr, minc:maxc])
        cropped_bin.append(mask[minr:maxr, minc:maxc])

    return zip(cropped_img, cropped_bin)


def grouper(iterable, n, fillvalue=None):
    '''Collect data into fixed-length chunks or blocks
    grouper('ABCDEFG', 3, 'x') --> ABC DEF Gxx'''
    # Adapted from: https://docs.python.org/3/library/itertools.html
    args = [iter(iterable)] * n
    return zip_longest(*args, fillvalue=fillvalue)

def remove_outliers(results):
    clean_results = {}
    outliers = {}
    for key ,val in results.items():
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

    def G(Z): return np.abs(Z).argmax()  # returns indices of min values
    def t_crit(N): return stats.t.isf(alpha / (2. * N), N - 2)
    def G_crit(N): return (N - 1.) / np.sqrt(N) * \
        np.sqrt(t_crit(N)**2 / (N - 2 + t_crit(N)**2))

    outliers = np.array([])
    while np.amax(np.abs(Z)) > G_crit(N):
        # update the outliers
        outliers = np.r_[outliers, X[G(Z)]]
        # remove outlier from array
        X = np.delete(X, G(Z))
        # repeat Z score
        Z = stats.zscore(X, ddof=1)
        N = len(X)

    return X, outliers

if __name__ == '__main__':
    main(sys.argv)
