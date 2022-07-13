#!/usr/bin/env python3

import time
import argparse
import npmcv


def run():
    # Argparse
    parser = argparse.ArgumentParser(prog='npmcv')
    parser.add_argument('-D', '--no_imgs', action='store_false',
                        help='Don\'t save overlays and individual cell images')
    parser.add_argument('-n', '--name', default=None, type=str,
                        nargs='?', help='Set basename for the output files')
    parser.add_argument('-V', '--version', action='version',
                        version='npmcv version: {}'.format(npmcv.__version__))
    parser.add_argument('path', type=str, metavar='<path>',
                        help='directory containing lif images')

    args = parser.parse_args()

    t = time.process_time()
    npmcv.main(vars(args))
    elapsed_time = (time.process_time() - t) / 60
    print('\x1b[1;32m' +
          'Total Runtime: {:.4} mins.'.format(elapsed_time) + '\x1b[0m \n')


if __name__ == '__main__':
    run()
