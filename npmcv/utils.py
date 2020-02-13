import re
from operator import itemgetter
from itertools import zip_longest, groupby


'''
Some functions that are currently not used
but might be useful
'''

def dvextractor(directory):
    '''
    directly extracts images from the deltavision file type
    requires mirc
    '''
    import mrc
    import numpy as np

    img_files = sorted([f for f in os.listdir(directory)
                        if f.endswith('.dv')], key=lambda f: f.lower())

    for im in img_files:
        dv = mrc.imread(im) 

        # Expecting a 4D numpy array (z,c,x,y) 
        # the second number indicates the "channel" 0 = dapi, 1 = npm  

        raw_dapi = im[0][0]
        raw_npm1 = im[0][1]


def img_counter(imgs):
    '''Determines how many Images were taken per slide
    to later combine those cells into one column

    Returns list containing slicing indexes to apply to raw data
    '''
    seq = []
    for img in imgs:
        m = re.search(r'.*Image(\d{3})_s00', img)
        if m:
            seq.append(int(m.group(1)))

    count = [list(map(itemgetter(1), g))
             for k, g in groupby(enumerate(seq), lambda x: x[0] - x[1])]
    ips = []
    for n in count:
        ips.append(len(n))

    start, stop = 0, ips[0]
    results = [[start, stop]]
    for n in ips[1:]:
        start = stop
        stop = start + n
        results.append([start, stop])

    return results


def sorter(raw_data, slice_index):
    '''Processing Raw Data'''

    print('Formating Data...\n')
    df = {}
    for stsp in slice_index:
        lst = []
        for key, val in list(raw_data.items())[stsp[0]:stsp[1]]:
            lst.extend(val)
        clst = [x for x in lst if not np.isnan(x)]
        slide = re.split(r'([\s]+)', key)[0]
        df[slide] = clst

    print('\x1b[6;30;42m' +
          'Done!: {n} Images Analyzed'.format(n=len(raw_data)) + '\x1b[0m', end='')

    return df


def verbosity(img, labeled, dist, relabeled):
    '''
    Helper function for graphing segmentation example
    '''
    fig, axes = plt.subplots(2, 2, figsize=(8, 8))
    ax = axes.flatten()
    titles = ["threshold_li", "distance", "watershed", "relabeled"]
    images = [img, labeled, dist, relabeled]
    for x, img, name in zip(range(4), images, titles):
        ax[x].imshow(img, cmap=plt.cm.nipy_spectral)
        ax[x].set_title(name)
        ax[x].set_axis_off()

    plt.tight_layout()
    plt.show()
    #plt.savefig('test.png', format='png')
