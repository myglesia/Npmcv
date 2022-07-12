import sys
import os
from itertools import zip_longest
from collections import defaultdict

import mahotas as mh
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from scipy import stats
from scipy import ndimage as ndi
from skimage import exposure, color
from skimage.util import img_as_float
from skimage.measure import regionprops
from skimage.filters import gaussian, threshold_li
from readlif.reader import LifFile



def main(argv):
    '''npmcv <directory>

    where the directory contains a Lecia .lif images
    This is also where the data is saved.
    '''
    directory = argv['path']
    os.chdir(directory)

    lif_files = [os.path.join(directory, f) for f in os.listdir(directory) if f.endswith(".lif") and not f.startswith('.')]

    if not lif_files: print("No LIF files found!")

    filecv = defaultdict(list)
    invcv = defaultdict(list)
    for lif in lif_files:

        imgcv = liffr(lif) # [('img', [CV]), ('img', [CV]), ... ]

        for k,v in imgcv:
            filecv[os.path.basename(lif)].extend(v)
            # {'file': [CVs], ...}
            invcv[k].extend(v)
            # {'img': [CVs], ...}

    results, outliers = remove_outliers(filecv)
    # handling save files
    # the name of the input folder
    if argv['output']:
        exp_name = os.path.basename(argv['output'])
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

def liffr(file):
    '''Processes a single .lif image.'''

    new = LifFile(file)
    fname = os.path.basename(new.filename)
    print('\x1b[4m' + 'Processing: {} '.format(fname) + '\x1b[0m')

    # Create a list of images using a generator
    img_list = [i for i in new.get_iter_image()]
    results = []
    for img in img_list:
        iname = os.path.splitext(fname)[0] + " " + img.name
        print(' {} '.format(img.name))
        try: # Check if image contains the correct number of channels
            # Returns Pillow objects
            ch_list = [i for i in img.get_iter_c(t=0, z=0)]
            dapi = ch_list[0]
            npm1 = ch_list[1]

        except IndexError:
            print('Incorrect number of channels in {}, skipping image...'.format(iname))
            continue

        img_cv = sip(np.array(dapi), np.array(npm1), iname)

        # each image is a column of cv
        # [('img', [CV]), ('img', [CV]), ... ]
        results.append((iname, img_cv))

    #raw_results
    return results

def sip(raw_dapi, raw_npm1, name):
    '''Single image Processing - The Images are actually opened here

    Returns
    -------
        results:    list
            list of stats for one image
    '''
    label = cell_segmentation(raw_dapi)

    results = []
    # no usuable cells
    if label is None:
        results.append(np.nan)
    else:
        raw_npm1 = exposure.rescale_intensity(img_as_float(raw_npm1))
        check(label, raw_npm1, name)
        cells = img2cells(label, raw_npm1)
        for i, (cell, mask) in enumerate(cells):
            cv = stats.variation(cell[mask], axis=None)
            verb(i, cell, mask, name)
            results.append(cv)
    return results


def cell_segmentation(img):
    '''Prepares a mask labeling each cell in the DAPI images,
    removing artifacts and cells touching the edges.

    Returns
    -------
    cleaned:    label
        labeled binary image for a single slide
    '''
    # Li Threshold...
    im = gaussian(img_as_float(img), sigma=1)  # time!
    bin_image = (im > threshold_li(im))
    labeled, nr_objects = mh.label(bin_image)

    # NOTE Change these values first
    cleaned, n_left = mh.labeled.filter_labeled(
        labeled, remove_bordering=True, min_size=10000, max_size=30000)

    print('\nNumber of cells counted: ' +
          '\x1b[3;31m' + '{}'.format(n_left) + '\x1b[0m')

    if n_left == 0:
        return None
    else:
        return cleaned


def verb(i, cell, mask, name):
    '''Verbose
    '''
    img = np.zeros(cell.shape, dtype=cell.dtype)
    img[:, :] = mask * cell

    os.makedirs('inv_cells', exist_ok=True)
    plt.imsave(os.path.join(os.getcwd(), 'inv_cells',
               '{0}_cell{1}.png'.format(name, i)), img)


def check(label, img, name):
    '''
    Quality check for dapi segmentation:
        Overlays the dapi segmentation onto full npm image
    '''

    rgbimg = color.label2rgb(label, img, alpha=0.2,
                             bg_label=0, bg_color=(0, 0, 0))

    os.makedirs('dapi_seg', exist_ok=True)
    newname = os.path.join(os.getcwd(), 'dapi_seg', '{0}_seg.png'.format(name))
    plt.imsave(os.path.join(os.getcwd(), 'dapi_seg',
               '{0}_seg.png'.format(name)), rgbimg)

    print('\nSegmentation Image saved: ' +
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


if __name__ == '__main__':
    main(sys.argv)
