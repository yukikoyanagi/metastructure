#!/usr/bin/env python3

import argparse
import pickle
import re
import pathlib
import numpy as np

def main(bp_file):
    with open(bp_file) as fh:
        data = fh.read()
    pattern = r'[0-9]{1,2}--[0-9]{1,2}:[PA]'
    ms = re.findall(pattern, data)
    pairs = []
    for m in ms:
        c = m.split(':')[1]
        f, s = m.split(':')[0].split('--')
        f = int(f)-1
        s = int(s)-1
        pairs.append((f,s,c))

    n = max([p[0] for p in pairs] + [p[1] for p in pairs]) + 1
    mat = np.zeros((n, n), dtype=int)
    for pair in pairs:
        i, j, c = pair
        if c == 'A':
            i, j = j, i
        if np.sum(mat[i,:] + mat[:,i] + mat[j,:] + mat[:,j]) > 1:
            continue
        else:
            mat[i,j] = 1

    p = pathlib.Path(bp_file)
    outf = p.with_suffix('.pkl')
    with open(outf, 'wb') as fh:
        pickle.dump(mat, fh)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('bp_file')
    args = parser.parse_args()
    main(args.bp_file)
