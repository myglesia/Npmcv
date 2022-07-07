import os
import numpy as np
from PIL import Image
from readlif.reader import LifFile


def _liffr(file):
    '''Processes a single .lif image.'''

    new = LifFile(file)
    fname = os.path.basename(new.filename)
    print('\x1b[4m' + 'Processing: {} '.format(fname) + '\x1b[0m')

    # Create a list of images using a generator
    img_list = [i for i in new.get_iter_image()]
    #results = {}
    results = []
    for img in img_list:
        print(' {} '.format(img.name))
        # Returns Pillow objects
        ch_list = [i for i in img.get_iter_c(t=0, z=0)]
        dapi = ch_list[0]
        npm1 = ch_list[1]

        iname = os.path.splitext(fname)[0] + " " + img.name
        img_cv = npmcv.sip(np.array(dapi), np.array(npm1), iname)

        # each image is a column of cv
        # [('img', [CV]), ('img', [CV]), ... ]
        results.append((iname, img_cv))

        #results[iname] = img_cv

    #raw_results
    return results