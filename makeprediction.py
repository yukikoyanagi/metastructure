#!/usr/bin/env python3
#

import argparse, json, logging, pickle
from itertools import permutations, product
import numpy as np

from fatgraph.fatgraph import FatgraphB

MAT_FILE = 'gbs_arr_20.pkl'

with open(MAT_FILE, 'rb') as fh:
    bg_mat = pickle.load(fh)

def partition(collection):
    if len(collection) == 1:
        yield [ collection ]
        return

    first = collection[0]
    for smaller in partition(collection[1:]):
        # insert `first` in each of the subpartition's subsets
        for n, subset in enumerate(smaller):
            yield smaller[:n] + [[ first ] + subset]  + smaller[n+1:]
        # put `first` in its own subset 
        yield [ [ first ] ] + smaller

def validorders(collection):
    "Return all valid reorderings of elements in collection. "
    "Collection is a list of lists. Rules for valid orders are; "
    "-within each list, the first element is smaller than the last"
    "-lists are ordered by the smallest element in each list"
    cs = sorted([sorted(c) for c in collection])
    ls = []
    for c in cs:
        l = []
        for p in permutations(c, len(c)):
            if p[0] < p[-1]:
                l.append(p)
        ls.append(l)
    return product(*ls)

def allpmats(n):
    "Yield all valid, complete pairing matrices of size n"
    for p in partition(list(range(n))):
        if min([len(s) for s in p]) < 2:
            continue
        for q in validorders(sorted(p)):
            #q is a tuple of tuples representing strand pairings
            links = [(i,j) for r in q for i,j in zip(r, r[1:])]
            for link_oris in product([True, False], repeat=len(links)):
                mat = np.zeros((n, n))
                for idx, par in zip(links, link_oris):
                    if par:
                        mat[idx] = 1
                    else:
                        mat[idx[::-1]] = 1
                yield mat

def main(pf, save, alpha, dbg):
    if dbg:
        logging.basicConfig(format='%(message)s', level=logging.DEBUG)
    else:
        logging.basicConfig(format='%(message)s', level=logging.INFO)

    with open(pf) as fh:
        pmat = np.array(json.load(fh))

    logging.debug(pmat)

    n = pmat.shape[0]
    hi_score = 0
    hi_mat = None
    beta = 1 - alpha
    for cmat in allpmats(n):
        if beta > 0:
            fg = FatgraphB.from_pmat(cmat)
            g = fg.genus
            b = len(fg.boundaries)
            k = max(len(v) for v in fg.vertices)//2
            try:
                topo_score = bg_mat[g,b,k] / np.sum(bg_mat[:,:,k])
            except IndexError as e:
                topo_score = 0
        else:
            topo_score = 0    
        score = alpha * np.sum(pmat*cmat) + beta * topo_score
        if score > hi_score:
            hi_score = score
            hi_mat = cmat

    logging.debug(hi_score)

    if save:
        with open(save, 'w') as fh:
            json.dump(hi_mat.tolist(), fh)
    else:
        print(hi_mat)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('pmat_file',
                        help='Pairing matrix json file.')
    parser.add_argument('-s', '--save',
                        help='Save output in this file.')
    parser.add_argument('-a', '--alpha', type=float, default=1.0,
                        help='Factor for betapro score. Beta is '
                        'computed as 1-alpha.')
    parser.add_argument('-d', '--debug', action='store_true',
                        help='Print debug messages.')
    args = parser.parse_args()
    main(args.pmat_file, args.save, args.alpha, args.debug)
