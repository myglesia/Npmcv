#!/usr/bin/env python

import sys
from npmcv import npmcv
__version__ = '0.3.0'

def run():
    for arg in sys.argv[1:]:
        npmcv.main(arg)

if __name__ == '__main__':
	run()
