#!/usr/bin/env python3

import argparse, pathlib, pickle, math
import numpy as np

def tririse(mat, orientation=True):
    "Make matrix upper-triangular by folding along the diagonal."
    "If orientation=True, multiply the lower-triangle by -1 before "
    "adding to the upper-triangle."
    upper = np.triu(mat, 1)
    lower = np.tril(mat, -1)
    if orientation:
        lower = -1 * lower
    return upper + lower.T

def compare(candidate, target, orientation=True, only=None, upto=None):
    "Return recall (TPR) and precision (PPV) of candidate "
    "matrix against target. If ignore_ori is True, ignore orientation."
    "If only (o) is given, only consider the o'th and -o'th diagonals."
    "If upto (u) is given, consider up to and including the u'th and "
    "-u'th diagonals."
    if only is not None:
        cmat = keeponly(candidate, only)
        tmat = keeponly(target, only)
    elif upto is not None:
        cmat = keepupto(candidate, upto)
        tmat = keepupto(target, upto)
    TP = np.sum(np.logical_and(cmat, tmat))
    FP = np.sum(np.where(cmat-tmat>0, 1, 0))
    FN = np.sum(np.where(tmat-cmat>0, 1, 0))
    if TP==0:
        TPR = 0
        PPV = 0
    else:
        TPR = TP/(TP+FN)
        PPV = TP/(TP+FP)
    return TPR, PPV

'''#Return accuracy
    cmat = tririse(candidate, orientation)
    tmat = tririse(target, orientation)
    res = np.where(cmat==tmat, 1, 0)
    n = res.shape[0]
    if only is not None:
        res = keeponly(res, only)
        m = n-only
    elif upto is not None:
        res = keepupto(res, upto)
        m = n*upto - upto*(upto+1)/2
    else:
        m = (n-1)*n/2
    return np.sum(np.triu(res, 1))/m
'''

def makeonematrix(pmat, omat):
    "Make a single pairing matrix from upper-triangular pairing "
    "matrix and orientation matrix. Anti-parallel connections are in "
    "the lower-triangular part. Negative entries are set to 0."
    pmat = np.triu(pmat, 1)
    pmat = np.where(pmat>0, pmat, 0)
    para_mat = np.where(omat, pmat, 0)
    anti_mat = np.where(np.logical_and(pmat, np.logical_not(omat)),
                        pmat, 0).T
    return para_mat + anti_mat

def keepupto(mat, u):
    "Keep up to u'th and -u'th diagonals (inclusive) of mat."
    return mat - np.triu(mat, u+1) - np.tril(mat, -u-1)

def keeponly(mat, o):
    "Keep only u'th and -u'th diagonals of mat"
    return np.diag(np.diag(mat, o), o) + np.diag(np.diag(mat, -o), -o)

def computeone(pf, tf, cutoff, only=None, upto=None):
    "Compute accuracy for a single candidate structure."
    pf = pathlib.Path(pf)
    with open(pf, 'rb') as fh:
        cand_pmat, cand_omat = pickle.load(fh)
    cand_mat = makeonematrix(cand_pmat, cand_omat)
    cand_mat = np.where(cand_mat<=cutoff, 0, 1)

    tf = pathlib.Path(tf)
    with open(tf, 'rb') as fh:
        mot = pickle.load(fh)
    targ_mat = np.zeros_like(cand_mat)
    for link in mot:
        i,j = sorted(link[0])
        if link[1][0] == link[1][1]:#'rr' or 'll'
            i,j = j,i
        targ_mat[i,j] = 1
    if upto is not None:
        cand_mat = keepupto(cand_mat, upto)
        targ_mat = keepupto(targ_mat, upto)
    elif only is not None:
        cand_mat = keeponly(cand_mat, only)
        targ_mat = keeponly(targ_mat, only)
    return compare(cand_mat, targ_mat, True, only, upto)

def main(pm, tdir, cutoff, only=None, upto=None):
    if pathlib.Path(pm).is_dir():
        acs = []
        for pf in pathlib.Path(pm).glob('*.pkl'):
            mf = pathlib.Path(tdir) / '{}.pkl'.format(pf.name.split('.')[0])
            acs.append(computeone(pf, mf, cutoff, only, upto))
        TPRs, PPVs = zip(*acs)
        print('Avg.TPR:{} Avg.PPV:{}'.format(
            sum(TPRs)/len(TPRs), sum(PPVs)/len(PPVs)))
    else:
        pf = pathlib.Path(pm)
        mf = pathlib.Path(tdir) / '{}.pkl'.format(pf.name.split('.')[0])
        TPR, PPV = computeone(pf, mf, cutoff, only, upto)
        print('TPR:{} PPV:{}'.format(TPR, PPV))

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Evaluation of alignment algorithm. If '
        'partial_matrix is a single file, returns reacall (TPR) and '
        'precision (PPV). If it is a directory, returns averages of '
        'recall and precision values.'
    )
    parser.add_argument('partial_matrix',
                        help='A pickle file containing partial '
                        'matrix, or a directory with such pickle files.')
    parser.add_argument('true_motif_dir',
                        help='Directory containing pickle files of '
                        'true motifs.')
    parser.add_argument('--cutoff', '-c', type=float, default=0,
                        help='Cutoff value for candidate pairing matrix.')
    group = parser.add_mutually_exclusive_group()
    group.add_argument('--only', '-o', type=int,
                       help="Only consider (upper and lower) o'th "
                       "diagonal.")
    group.add_argument('--upto', '-u', type=int,
                       help="Only consider up to and including "
                       "(upper and lower) u'th diagonal.")
    args = parser.parse_args()
    main(args.partial_matrix, args.true_motif_dir,
         args.cutoff, args.only, args.upto)
