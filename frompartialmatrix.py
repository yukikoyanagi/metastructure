#!/usr/bin/env python3

CUTOFF = 0.2

import argparse, pickle, logging, pathlib
from itertools import combinations, product, accumulate
from operator import mul
from datetime import datetime
import numpy as np

class BarrelError(Exception):
    pass

def combinations_repeat(it, r, n):
    "Return n-tuple of r-length subsequences of elements from input "
    "iterable."
    """
    combinations_repeat('ABCDE', 2, 2) --> (AB, CD) (AB, CE) (AB, DE)
                                           (AC, BD) (AC, BE) (AC, DE)
                                           (AD, BC) (AD, BE) (AD, CE)
                                           (AE, BC) (AE, BD) (AE, CD)
                                           (BC, DE) (BD, CE) (BE, CD)
    """
    indices = range(len(it))
    if r*n > len(it) or r*n==0:
        return
    for i in combinations(combinations(indices, r),n):
        #yield only if the result does not contain duplicates
        if len([k for j in i for k in j]) \
           == len({k for j in i for k in j}):
            res = []
            for j in i:
                subseq = tuple(it[k] for k in j)
                res.append(subseq)
            yield tuple(res)

def offdiagonal(s_mat, o_mat, n):
    motif = set()
    m = s_mat.shape[0]
    for i, j in zip(range(m), range(n, m)):
        if s_mat[i,j] < CUTOFF:
            continue
        else:
            motif.add(((i,j), o_mat[i,j]))
    return motif

def isbifurcated(mat):
    "Determine if given *full* pairing matrix has bifurcations."
    for row in mat:
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
    "Determine if given *full* pairing matrix has barrels."
    seen = []
    w_mat = np.where(mat>0, mat, 0)
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

def hasisolated(mat):
    "Check if given *full* pairing matrix results in an isolated strand"
    isolated = [i for i in range(len(mat))
                if np.count_nonzero((mat>0)[i])==0]
    return bool(isolated)

def makepartialmatrix(s_mat):
    "Partial pairing matrix from score matrix. "
    "Forbid barrels and bifurcations. Isolated strands are allowed."
    '''
    array([[ 0,.8,.5,-1, 0],       array([[ 0, 1,-1,-1, 0],
           [ 0, 0,.6,-1,-1],              [ 0, 0, 1,-1,-1],
           [ 0, 0, 0,.4,-1],   |-->       [ 0, 0, 0, 1,-1],
           [ 0, 0, 0, 0,.7],              [ 0, 0, 0, 0, 1],
           [ 0, 0, 0, 0, 0]])             [ 0, 0, 0, 0, 0]])
    '''
    w_mat = np.where(np.logical_and(s_mat>0, s_mat<CUTOFF), -1, s_mat)
    m_mat = np.zeros_like(w_mat, dtype=int)
    while w_mat.max() > 0:
        idx = np.unravel_index(np.argmax(w_mat), w_mat.shape)
        m_mat[idx] = 1
        r_idx = idx[::-1]
        m_mat[r_idx] = 1
        if isbifurcated(m_mat) or hasbarrel(m_mat):
            m_mat[idx] = 0
            m_mat[r_idx] = 0
        w_mat[idx] = 0
    mat = np.where(s_mat > 0, m_mat, s_mat)
    return np.where(np.logical_and(mat==0, s_mat!=0), -1, mat)

def completions(p_mat, o_mat):
    "All possible completions of partial pairing matrix, excluding "
    "bifurcations, barrels, and isolated strands."
    "Return a list of tuples consisting score and orientation matrices."
    w_mat = p_mat + p_mat.T
    #edges includes isolated strands
    edges = [i for i in range(len(w_mat))
             if np.count_nonzero((w_mat>0)[i])<2]
    #Isolated strands are 2-valent, so add them again
    isolated = [i for i in range(len(w_mat))
                if np.count_nonzero((w_mat>0)[i])==0]
    es = sorted(edges + isolated)
    if not any([isbifurcated(w_mat),
                hasbarrel(w_mat),
                hasisolated(w_mat)]):
        yield np.where(p_mat>0, p_mat, 0), o_mat
    for n in range(len(es)//2 + 1):
        for indices in set(combinations_repeat(es, 2, n)):
            if len(set(indices) & set(zip(*np.where(p_mat>0))))>0:
                #Can't pair existing pair
                continue
            if len([idx for idx in indices if idx[0]==idx[1]])>0:
                #An isolated strand is bonded to itself
                continue
            if -1 in [p_mat[idx] for idx in indices]:
                #Pairing is forbidden where p_mat is -1
                continue
            c_mat = np.copy(w_mat)
            for idx in indices:
                c_mat[idx] = 1
                r_idx = idx[::-1]
                c_mat[r_idx] = 1
            if any([isbifurcated(c_mat),
                    hasbarrel(c_mat),
                    hasisolated(c_mat)]):
                #c_mat results in invalid motif
                continue
            c_mat = np.where(c_mat>0, c_mat, 0)
            #c_mat is valid. Construct all possible omat
            x_mat = np.copy(o_mat)
            for rs in product([False, True], repeat=n):
                for i, v in enumerate(rs):
                    x_mat[indices[i]] = v
                yield np.triu(c_mat, 1), x_mat

def tomotif(p_mat, o_mat):
    "Make motif from valid pairing and orientation matrices"
    wmat = p_mat.copy()
    motif = set()
    while wmat.any():
        idx = np.unravel_index(np.argmax(wmat), wmat.shape)
        motif.add(((idx[0], idx[1]), o_mat[idx]))
        wmat[idx] = 0
    return motif

def main(partial, n, save, dbg):
    if dbg:
        logging.basicConfig(format='%(message)s', level=logging.DEBUG)
    else:
        logging.basicConfig(format='%(message)s', level=logging.INFO)
    logger = logging.getLogger()
    pf = pathlib.Path(partial)
    pid = pf.stem.split('.')[0].upper()
    start = datetime.now()
    logger.info('{}: Processing {}...'.format(start, pid))
    with open(partial, 'rb') as fh:
        s_mat, o_mat = pickle.load(fh)
    logger.debug(s_mat)
    logger.debug(o_mat)
    if n:
        motif = offdiagonal(s_mat, o_mat, n)
        if save:
            sd = pathlib.Path(save)
            of = sd / '{}.pkl'.format(pid)
            with open(of, 'wb') as fh:
                pickle.dump(motif, fh)
        else:
            print(motif)
    else:
        part_mat = makepartialmatrix(s_mat)
        logger.debug('Partial matrix:')
        logger.debug(part_mat)
        motifs = []
        for pmat, omat in completions(part_mat, o_mat):
            logger.debug(pmat)
            motifs.append(tomotif(pmat, omat))
        end = datetime.now()
        d = end - start
        logger.info('{}: Finished {} in {} seconds'.format(
            end, pid, d.seconds))
        if save:
            sd = pathlib.Path(save)
            of = sd / '{}.pkl'.format(pid)
            with open(of, 'wb') as fh:
                pickle.dump(motifs, fh)
        else:
            print(motifs)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Construct protein metastructure from partial '
        'motif by going through all possible completions of partial '
        'matrix and applying genus filter.')
    parser.add_argument('partial',
                        help='Partial motif pkl file.')
    parser.add_argument('-o', '--only', type=int,
                        help='Only consider o\'th diagonal, ignoring '
                        'validity of the resulting motif.')
    parser.add_argument('-s', '--save',
                        help='Save output to this directory.')
    parser.add_argument('-d', '--debug', action='store_true',
                        help='Print debug messages.')
    args = parser.parse_args()
    main(args.partial, args.only, args.save, args.debug)
