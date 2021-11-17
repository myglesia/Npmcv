#!/usr/bin/env python

import sys
import time
from npmcv import npmcv


def run():
    t = time.process_time()
    npmcv.main(sys.argv[1:])
    elapsed_time = (time.process_time() - t) / 60
    print('\x1b[1;32m' +
          'Total Runtime: {:.4} mins.'.format(elapsed_time) + '\x1b[0m \n')


if __name__ == '__main__':
    run()
