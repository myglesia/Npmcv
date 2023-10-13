import os
from itertools import zip_longest
from collections import defaultdict

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from scipy import stats
from scipy import ndimage as ndi

from skimage import exposure, color
from skimage.util import img_as_float
from skimage.segmentation import clear_border, watershed
from skimage.measure import label, regionprops
from skimage.filters import gaussian, threshold_otsu
from skimage.morphology import opening, disk
from readlif.reader import LifFile


def main(argv):
    '''npmcv <directory>

    where the directory contains a Lecia .lif images
    This is also where the data is saved.
    '''
    directory = os.path.abspath(argv['path'])

    os.chdir(directory)

    lif_files = [os.path.join(directory, f) for f in os.listdir(
        directory) if f.endswith(".lif") and not f.startswith('.')]

    if not lif_files:
        print("No LIF files found!")

    filecv = defaultdict(list)
    invcv = defaultdict(list)

    # [('img', [CVs]), ('img', [CVs]), ... ]
    for lif in lif_files:
        imgcv = liffr(lif, argv) 
        for k, v in imgcv:
            # {'file': [CVs], ...}
            filecv[os.path.basename(lif)].extend(v)
            # {'img': [CVs], ...}
            invcv[k].extend(v)

    results, outliers = remove_outliers(filecv)
    
    # handling save files
    if argv['name']:
        exp_name = os.path.basename(argv['name'])
    else:
        exp_name = os.path.basename(os.path.abspath(directory))

    def dtdf(d): return pd.DataFrame(
        dict([(k, pd.Series(v)) for k, v in d.items()]))

    df = dtdf(results)
    dtdf(results).to_csv('{}.csv'.format(exp_name), index=False)
    dtdf(invcv).to_csv('{}_RAW.csv'.format(exp_name), index=False)
    dtdf(outliers).to_csv('{}_OUT.csv'.format(exp_name), index=False)

    # Final Stats
    print('\n')
    print('\x1b[1;36m' +
          'Experiment {0} Completed!.'.format(exp_name) + '\x1b[0m')
    print('\x1b[1;36m' + ' {n} '.format(n=str(sum(df.count()))) + '\x1b[0m' +
          'Total Cells Counted.')


def liffr(file, argv):
    '''Processes a single .lif image.'''

    save=argv['no_imgs']

    lif_img = LifFile(file)
    fname = os.path.basename(lif_img.filename)
    print('\n\x1b[4m' + 'Working on... {} '.format(fname) + '\x1b[0m')

    # Create a list of images in Lif file using a generator
    img_list = [i for i in lif_img.get_iter_image()]
    results = []
    for img in img_list:
        iname = os.path.splitext(fname)[0] + " " + img.name
        print('\n    Processing: {} '.format(img.name))
        if img.bit_depth[0] < 16:
            print('\x1b[1;31m' + '    Warning {}-bit image.'.format(img.bit_depth[0]) + '\x1b[0m')
        try:  
            # Check if image contains the correct number of channels
            ch_list = [i for i in img.get_iter_c(t=0, z=0)] # Returns Pillow objects
            dapi = ch_list[0]
            npm1 = ch_list[1]

        except IndexError:
            print('incorrect number of channels, skipping {}...'.format(iname))
            continue

        img_cv = sip(np.array(dapi), np.array(npm1), iname, argv)

        # each image is a column of cv
        # [('img', [CV]), ('img', [CV]), ... ]
        results.append((iname, img_cv))

    # Raw_results
    return results


def sip(raw_dapi, raw_npm1, name, argv):
    '''Single image Processing - The Images are actually processed here


    Returns
    -------
        results:    list
            list of stats for one image
    '''

    # Prepares a mask labeling each cell in the DAPI images, after
    # removing artifacts and cells touching the edges.

    smooth = gaussian(raw_dapi, sigma=1.2)
    smooth = opening(smooth, disk(20))

    thresh = (smooth > threshold_otsu(smooth))
    fill = ndi.binary_fill_holes(thresh)

    # Identify some background and foreground pixels from the intensity values.
    # These pixels are used as seeds for watershed.
    markers = np.zeros_like(smooth)
    foreground, background = 1, 2
    markers[fill == True] = foreground
    markers[fill == False] = background

    ws = watershed(smooth, markers)

    # remove artifacts connected to image border and label image regions
    label_image = label(clear_border(ws) == foreground)

    # remove labels that are too large or too small (based on pixel area)
    d_label, n_left = clear_size(label_image, argv)

    print('    Number of cells counted: ' +
          '\x1b[3;32m' + '{}'.format(n_left) + '\x1b[0m')

    save=argv['no_imgs']
    results = []
    if n_left == 0:
        # no usuable cells
        results.append(np.nan)
    else:
        #raw_npm1 = exposure.rescale_intensity(img_as_float(raw_npm1))
        if save:
            check(d_label, raw_npm1, name)
        cells = img2cells(d_label, raw_npm1)
        for i, (cell, mask) in enumerate(cells):
            cv = stats.variation(cell[mask], axis=None)
            if save:
                verb(i, cell, mask, name)
            results.append(cv)
    return results


def verb(i, cell, mask, name):
    '''Verbose
    '''
    img = np.zeros(cell.shape, dtype=cell.dtype)
    img[:, :] = mask * cell

    os.makedirs('inv_cells', exist_ok=True)
    plt.imsave(os.path.join(os.getcwd(), 'inv_cells',
                            '{0}_cell{1}.png'.format(name, i)), img)


def check(label, img, name):
    '''Quality check for dapi segmentation:
        Overlays the dapi segmentation onto full npm image
    '''

    rgbimg = color.label2rgb(label, img, alpha=0.2,
                             bg_label=0, bg_color=(0, 0, 0))

    os.makedirs('dapi_seg', exist_ok=True)
    newname = os.path.join(os.getcwd(), 'dapi_seg', '{0}_seg.png'.format(name))
    plt.imsave(os.path.join(os.getcwd(), 'dapi_seg',
                            '{0}_seg.png'.format(name)), rgbimg)

    print('    Segmentation Image saved: ' +
          '\x1b[4m' + '{}\n'.format(newname) + '\x1b[0m')


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

def clear_size(labeled_img, argv):
    '''Removes labels outside of the min/max area size. 

    Returns
    -------
        out:   relabeled image, array-like  
        nleft:  number of remaining labels (cells)
    '''
    min_size=argv['min_size']
    max_size=argv['max_size']
    # Re-label, in case we are dealing with a binary out
    nlabels, number = label(labeled_img, background=0, return_num=True)
    
    borders_indices = [] # list of labels outside area range.
    for region in regionprops(nlabels):
        if region.area >= max_size or region.area < min_size:
            borders_indices.append(region.label)
    indices = np.arange(number + 1)  # 0 (background) is not counted as a label in regionprop

    # returns an array 'label_mask' as True at position of labels outside area range
    label_mask = np.isin(indices, np.unique(borders_indices))
    # create mask for pixels to clear
    mask = label_mask[nlabels.reshape(-1)].reshape(nlabels.shape)
    nlabels[mask] = 0

    #re-label image
    out, nleft = label(nlabels, background=0, return_num=True)
    return out, nleft

# Utility and stats functions
def grouper(iterable, n, fillvalue=None):
    '''Collect data into fixed-length chunks or blocks
    grouper('ABCDEFG', 3, 'x') --> ABC DEF Gxx

    Adapted from: https://docs.python.org/3/library/itertools.html
    '''
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
