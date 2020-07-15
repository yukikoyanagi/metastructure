#!/usr/bin/env python3

import argparse, pickle, pathlib
import numpy as np

def main(f, t, save):
    with open(f, 'rb') as fh:
        data = pickle.load(fh)
    if t=='m':
        smat, omat = data
        #print(smat, omat)
        pmat = np.zeros_like(smat)
        for i, j in np.ndindex(pmat.shape):
            if i<j:
                if omat[i,j]:
                    pmat[i,j] = smat[i,j]
                    pmat[j,i] = -1
                else:
                    pmat[j,i] = smat[i,j]
                    pmat[i,j] = -1
        #print(pmat)
        pmat = np.triu(np.tril(pmat, 5), -5)
        #print(pmat)

    if save:
        fn = pathlib.Path(f)
        of = pathlib.Path(save) / fn.with_suffix('.pkl').name
        with open(of, 'wb') as fh:
            pickle.dump(pmat, fh)
    else:
        print(pmat)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('input_file')
    parser.add_argument('-t', '--input_type', default='m',
                        help='Type of input file. Pairing & orientation '
                        'matrix pair: m.')
    parser.add_argument('-s', '--save',
                        help='Save output in this directory.')
    args = parser.parse_args()
    main(args.input_file, args.input_type, args.save)
