#!/usr/bin/env python3
#
# Produce strand pairing matrix from betapro output (from
# predict_beta_ss.sh). Written because predict_strand.sh does not
# produce full pairing matrix, but a symmetric one.

import argparse, re, json, logging
from itertools import combinations
import numpy as np

def getscore(mat, s1, s2, parallel):
    if len(s1) < len(s2):
        s = s1
        t = s2
    else:
        s = s2
        t = s1
    if not parallel:
        t = t[::-1]
    score = 0
    for k in range(len(t) - len(s) + 1):
        newscore = sum(mat[i,j] for i,j in zip(s, t[k:]))
        if newscore > score:
            score = newscore
    return score

def main(bpf, save, bridge, dbg):
    if dbg:
        logging.basicConfig(format='%(message)s', level=logging.DEBUG)
    else:
        logging.basicConfig(format='%(message)s', level=logging.INFO)

    strands = []
    rawm = []
    with open(bpf) as fh:
        for idx, line in enumerate(fh):
            if idx == 2:
                raw_ss = line.strip()
            if idx > 4:
                rawm.append(line.strip())
    rawm = ' '.join(rawm)
    logging.debug(raw_ss)
    logging.debug(rawm)

    ms = re.finditer(r'E+', raw_ss)
    n = 0
    for m in ms:
        std = list(range(n, n+m.end()-m.start()))
        strands.append(std)
        n += len(std)
    n = len([i for std in strands for i in std])
    if not bridge:
        strands = [s for s in strands if len(s)>1]
    logging.debug(strands)


    mat = np.fromstring(rawm, sep=' ').reshape((n,n))
    logging.debug(mat)

    m = len(strands)
    pmat = np.zeros((m, m))
    for i, j in combinations(range(m), 2):
        pmat[i,j] = getscore(mat, strands[i], strands[j], True)
        pmat[j,i] = getscore(mat, strands[i], strands[j], False)
    if save:
        with open(save, 'w') as fh:
            json.dump(pmat.tolist(), fh)
    else:
        print(pmat)

if __name__=='__main__':
    parser = argparse.ArgumentParser(
        description='Produce pairing matrix from BetaPro '
        'predict_beta_ss.sh output.'
    )
    parser.add_argument('bp_file',
                        help='BetaPro output containing residue '
                        'pairing matrix (but not strand pairings).')
    parser.add_argument('-s', '--save',
                        help='Save output in this json file.')
    parser.add_argument('-b', '--bridge', action='store_true',
                        default=False,
                        help='Count single beta residue as strand.')
    parser.add_argument('-d', '--debug', action='store_true',
                        help='Print debug messages.')
    args = parser.parse_args()
    main(args.bp_file, args.save, args.bridge, args.debug)
