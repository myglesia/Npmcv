#!/usr/bin/env python

import sys
from npmcv import npmcv
__version__ = '0.3.0'

def run():
	args = sys.argv
	npmcv.main(args)

if __name__ == '__main__':
	run()
