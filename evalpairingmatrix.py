#!/usr/bin/env python3

import argparse, pathlib, pickle, math, logging
import numpy as np

class BarrelError(Exception):
    pass

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
    else:
        cmat = candidate
        tmat = target
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

def keepupto(mat, u):
    "Keep up to u'th and -u'th diagonals (inclusive) of mat."
    return mat - np.triu(mat, u+1) - np.tril(mat, -u-1)

def keeponly(mat, o):
    "Keep only u'th and -u'th diagonals of mat"
    return np.diag(np.diag(mat, o), o) + np.diag(np.diag(mat, -o), -o)

def computeone(cand_mat, targ_mat, only=None, upto=None):
    "Compute accuracy for a single candidate structure."
    if upto is not None:
        cand_mat = keepupto(cand_mat, upto)
        targ_mat = keepupto(targ_mat, upto)
    elif only is not None:
        cand_mat = keeponly(cand_mat, only)
        targ_mat = keeponly(targ_mat, only)
    return compare(cand_mat, targ_mat, True, only, upto)

def isbifurcated(mat):
    "Determine if given pairing matrix has bifurcations."
    wmat = mat + mat.T
    for row in wmat:
        if np.count_nonzero(row>0) > 2:
            return True
    return False

def findsheetat(k, mat):
    if not np.allclose(mat, mat.T):
        raise ValueError('Not a symmetric matrix.')
    sheet = [k]
    row = mat[k]
    if np.count_nonzero(row) > 1:
        raise ValueError
    elif np.count_nonzero(row) < 1:
        return sheet
    prevk = k
    thisk = np.flatnonzero(row)[0]
    thisrow = mat[thisk]
    while np.count_nonzero(thisrow) != 1:
        if thisk in sheet:
            raise BarrelError
        sheet.append(thisk)
        i, j = np.flatnonzero(thisrow)
        if i == prevk:
            nextk = j
        else:
            nextk = i
        prevk, thisk = thisk, nextk
        thisrow = mat[thisk]
    else:
        if thisk in sheet:
            raise BarrelError
        sheet.append(thisk)
    return sheet

def hasbarrel(mat):
    "Determine if given pairing matrix has barrels."
    wmat = mat + mat.T
    seen = []
    w_mat = np.where(wmat>0, wmat, 0)
    for i, row in enumerate(w_mat):
        if i not in seen and np.count_nonzero(row) < 2:
            try:
                sheet = findsheetat(i, w_mat)
            except BarrelError as e:
                return True
            seen.extend(sheet)
    if len(seen) < len(w_mat):
        return True
    return False

def makepmat(smat):
    "Make pairing matrix (0 or 1) from score matrix"
    wmat = np.copy(smat)
    pmat = np.zeros_like(wmat, dtype=int)
    if len(pmat) == 2:
        if wmat[0,1] > wmat[1,0]:
            pmat[0,1] = 1
        else:
            pmat[1,0] = 1
        return pmat
    else:
        while wmat.max() > 0:
            idx = np.unravel_index(np.argmax(wmat), wmat.shape)
            r_idx = idx[::-1]
            pmat[idx] = 1
            if isbifurcated(pmat) or hasbarrel(pmat):
                pmat[idx] = 0
            wmat[idx] = 0
            wmat[r_idx] = 0
        return pmat

def main(pm, tdir, only=None, upto=None, dbg=False):
    if dbg:
        logging.basicConfig(format='%(message)s', level=logging.DEBUG)
    else:
        logging.basicConfig(format='%(message)s', level=logging.INFO)
    logger = logging.getLogger()

    if pathlib.Path(pm).is_dir():
        acs = []
        for pf in pathlib.Path(pm).glob('*.pkl'):
            logger.debug('Processing {}...'.format(pf.name))
            with open(pf, 'rb') as fh:
                mat = pickle.load(fh)
            pmat = makepmat(mat)
            mf = pathlib.Path(tdir) / '{}.pmat.pkl'.format(pf.name.split('.')[0])
            with open(mf, 'rb') as fh:
                mmat = pickle.load(fh)
            acs.append(computeone(pmat, mmat, only, upto))
        TPRs, PPVs = zip(*acs)
        print('Avg.TPR:{} Avg.PPV:{}'.format(
            sum(TPRs)/len(TPRs), sum(PPVs)/len(PPVs)))
    else:
        pf = pathlib.Path(pm)
        logger.debug('Processing {}...'.format(pf.name))
        with open(pf, 'rb') as fh:
            mat = pickle.load(fh)
        logger.debug(mat)
        pmat = makepmat(mat)
        logger.debug(pmat)
        mf = pathlib.Path(tdir) / '{}.pmat.pkl'.format(pf.name.split('.')[0])
        with open(mf, 'rb') as fh:
            mmat = pickle.load(fh)
        logger.debug(mmat)
        TPR, PPV = computeone(pmat, mmat, only, upto)
        print('TPR:{} PPV:{}'.format(TPR, PPV))

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Evaluation of candidate pairing matrix. If '
        'cand_matrix is a single file, returns reacall (TPR) and '
        'precision (PPV). If it is a directory, returns averages of '
        'recall and precision values.'
    )
    parser.add_argument('cand_matrix',
                        help='A pickle file containing candidate '
                        'matrix, or a directory with such pickle files.')
    parser.add_argument('true_motif_dir',
                        help='Directory containing pickle files of '
                        'true pairing matrices.')
    group = parser.add_mutually_exclusive_group()
    group.add_argument('--only', '-o', type=int,
                       help="Only consider (upper and lower) o'th "
                       "diagonal.")
    group.add_argument('--upto', '-u', type=int,
                       help="Only consider up to and including "
                       "(upper and lower) u'th diagonal.")
    parser.add_argument('--debug', '-d', action='store_true',
                        help='Print debug messages.')
    args = parser.parse_args()
    main(args.cand_matrix, args.true_motif_dir,
         args.only, args.upto, args.debug)
