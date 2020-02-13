#!/usr/bin/env python3

import sys
import os
import re
import mrc
from itertools import zip_longest
from collections import defaultdict

import mahotas as mh
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from PIL import Image
from scipy import stats, misc
from scipy import ndimage as ndi
from skimage.util import img_as_float
from skimage.measure import regionprops
from skimage.filters import gaussian, threshold_li

from skimage import color

def dv(img_files):
    '''
    return raw results: {'img': [CVs]}
    '''

    results = {}

    for im in img_files:
        dv_img = mrc.imread(im)
    # Expecting a 4D numpy array (z,c,x,y)
    # the second number indicates the "channel" 0 = dapi, 1 = npm
        raw_dapi = dv_img[0][0]
        raw_npm1 = dv_img[0][1]

        img_cv = sip(raw_dapi, raw_npm1, os.path.basename(im))
        # single image per column
        results[os.path.basename(im)] = img_cv

    return results

    # This regex stuff can't be done unless original file names are strictly controlled
    #pattern = re.compile(r"(?P<trial>^.+)\.(?P<ext>.+$)")

def tiff(img_files):
    '''
    return raw results: {'img': [CVs]}
    '''

    results = {}

    for dapi, npm1 in grouper(img_files, 2):
        print('\x1b[4m' + 'Processing: {}'.format(npm1) + '\x1b[0m')
        with Image.open(dapi) as d, Image.open(npm1) as n:
            raw_dapi = np.array(Image.open(dapi))
            raw_npm1 = np.array(Image.open(npm1))

            img_voc = sip(raw_dapi, raw_npm1, os.path.basename(im))
            # single image per column
            results[os.path.basename(im)] = img_cv

    return results

def main(argv):
    '''npmcv <directory>
    DAPI and NPM1 channels from a microscope image must be split into separate tif files. Files should be listed sequentially.

    e.g. "Image01.dapi.tif", "Image01.npm.tif", "Image02.dapi.tif",
    "Image02.npm.tif",...
    '''
    directory = argv[0]
    #'.dv'  or '.tiff'
    filetype = argv[1]
    # move cwd to where the images are
    # this is where everything will save
    os.chdir(directory)    

    # file list sorted alphabetically, case-insensitive
    # if .tiff
    # Like the form:
    #   ActD.lif - Image002_ch01.tif
    #   ActD.lif - Image002_ch02.tif
    img_files = sorted([os.path.join(directory, f) for f in os.listdir(directory)
                        if f.endswith(filetype) and not f.startswith('.')], key=lambda f: f.lower())

    raw_results = dv(img_files)

    # BETTER MATCHING
    #pattern = re.compile(r"(?P<slide>.*)\.*")
    #raw_results = defaultdict(list)
    #r = pattern.match(npm1)
    #raw_results[r.group(1)].extend(img_voc)

    results, outliers = remove_outliers(raw_results)

    # handling save files
    exp_name = os.path.basename(os.path.abspath(directory))

    def dtdf(d): return pd.DataFrame(
        dict([(k, pd.Series(v)) for k, v in d.items()]))

    df = dtdf(results)
    dtdf(results).to_csv('{}.csv'.format(exp_name), index=False)
    dtdf(raw_results).to_csv('{}_RAW.csv'.format(exp_name), index=False)
    dtdf(outliers).to_csv('{}_OUT.csv'.format(exp_name), index=False)

    # Final Stats
    print('\n')
    print('\x1b[1;36m' +
          ' Experiment {0} Completed!.'.format(exp_name) + '\x1b[0m')
    #print('\x1b[1;30m' +' {n} '.format(n=len(sep_results)) + '\x1b[0m' + 'Images Analyzed,', end='')
    print('\x1b[1;36m' + ' {n} '.format(n=str(sum(df.count()))) + '\x1b[0m' +
          'Cells Counted.\n')


def sip(raw_dapi, raw_npm1, name):
    '''Single image Processing - The Images are actually opened here

    Returns
    -------
        results:    list
            list of stats for one image
    '''
        # mrc deltavision function can insert here
    label = cell_segmentation(raw_dapi)

    results = []
    # no usuable cells
    if label is None:
        results.append(np.nan)
    else:
        check(label, raw_npm1, name)
        cells = img2cells(label, raw_npm1)
        print("Calculating CVs...")
        for i, (cell, mask) in enumerate(cells):
            cv = stats.variation(cell[mask], axis=None)
            verb(i, cell, mask, name)
            results.append(cv)

    return results

def cell_segmentation(img):
    '''Prepares a mask labeling each cell in the DAPI images,
    removing artifacts and cut off cells.

    Returns
    -------
    cleaned:    label
        labeled binary image for a single slide
    '''

    print("Segmentation Calculation...")
    # Li Threshold...
    im = gaussian(img_as_float(img), sigma=1)  # time!
    bin_image = (im > threshold_li(im))
    labeled, nr_objects = mh.label(bin_image)

    # NOTE Change these values first
    cleaned, n_left = mh.labeled.filter_labeled(
        labeled, remove_bordering=True, min_size=10000, max_size=30000)

    print('Final Number of cells: ' +
          '\x1b[3;31m' + '{}\n'.format(n_left) + '\x1b[0m')

    if n_left == 0:
        return None
    else:
        return cleaned


def verb(i, cell, mask, name):
    img = np.zeros(cell.shape, dtype=cell.dtype)
    img[:, :] = mask * cell

    #name = os.path.dirname(ppath) + '/in_cells/{0}_c{1}.png'.format(os.path.basename(ppath), i)
    #os.makedirs(os.path.dirname(name), exist_ok=True)

    os.makedirs('inv_cells', exist_ok=True)
    plt.imsave(os.path.join(os.getcwd(),'inv_cells','{0}_cell{1}.png'.format(name, i)), img)

def check(label, img, name):
    '''
    Quality check for dapi segmentation:
        Overlays the dapi segmentation onto full npm image
    '''

    rgbimg = color.label2rgb(label, img, alpha=0.2,
                             bg_label=0, bg_color=(0, 0, 0))


    os.makedirs('dapi_seg', exist_ok=True)
    plt.imsave(os.path.join(os.getcwd(),'dapi_seg','{0}seg.png'.format(name)), rgbimg)

    print('Segmentation Image saved: ' +
          '\x1b[4m' + '{}\n'.format(name) + '\x1b[0m')


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
    for key, val in results.items():
        val = [x for x in val if not np.isnan(x)]  # remove 'NaN'
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


if __name__ == '__main__':
    main(sys.argv)
