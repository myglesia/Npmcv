#!/usr/bin/env python

import sys
import time
from npmcv import npmcv


def run():
    t = time.process_time()
#    for arg in sys.argv[1:]:
#        npmcv.main(arg)
    npmcv.main(sys.argv[1:])
    elapsed_time = (time.process_time() - t) / 60
    print('\x1b[1;32m' +
          'Total Runtime: {:.4} mins.'.format(elapsed_time) + '\x1b[0m')


if __name__ == '__main__':
    run()
