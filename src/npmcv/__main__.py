#!/usr/bin/env python3

import time
import argparse
import Npmcv

def run():
    # Argparse
    parser = argparse.ArgumentParser(prog='npmcv')
    parser.add_argument('-o', '--output', default=None, type=str, nargs='?', help='output file name',)
    parser.add_argument('-V', '--version', action='version', version='npmcv version: {}'.format(Npmcv.__version__))
    parser.add_argument('path', type=str, help='directory containing lif images', metavar='<path>')
    args = parser.parse_args()


    t = time.process_time()
    Npmcv.main(vars(args))
    elapsed_time = (time.process_time() - t) / 60
    print('\x1b[1;32m' +
          'Total Runtime: {:.4} mins.'.format(elapsed_time) + '\x1b[0m \n')

if __name__ == '__main__':
    run()