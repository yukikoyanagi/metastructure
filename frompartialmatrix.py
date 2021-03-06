#!/usr/bin/env python3

MAT_FILE = 'gbs_arr.pkl'

import argparse, pickle, logging, pathlib
from itertools import combinations, product, accumulate
from operator import mul
from datetime import datetime
import numpy as np

from fatgraph.fatgraph import Fatgraph

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

def findsheets(p_mat):
    "Sheets from pairting matrix"
    '''
    findsheets([[0,0,1,0],
                [0,0,1,0],
                [0,0,0,0],
                [0,0,0,0]])
    --> {(0,2,1),(3,)}
    '''
    wmat = np.triu(p_mat,1) + np.triu(p_mat,1).T
    edges = [i for i in range(len(wmat))
             if np.count_nonzero(wmat[i])<2]
    sheets = set()
    while edges:
        e = edges.pop(0)
        sheet = tuple(findsheetat(e, wmat))
        if sheet[0] > sheet[-1]:
            sheet = sheet[::-1]
        sheets.add(sheet)
        if sheet[-1] in edges:
            edges.remove(sheet[-1])
    return sheets

def makevertices(sheets, o_mat):
    "Vertices from sheets and orientation matrix."
    '''
    [[0,1,2]], [[False, False, False],
                [False, False, True],
                [False, False, False]]
    --> [[(0,1),(11,10),(21,20)]]
    '''
    omat = np.triu(o_mat, 1) + np.triu(o_mat, 1).T
    vertices = []
    for sheet in sheets:
        oseq = []
        for i, j in zip(sheet, sheet[1:]):
            if omat[i,j]:
                oseq.append(1)
            else:
                oseq.append(-1)
            oseq = list(accumulate(oseq, mul))
        vertex = []
        i = sheet[0]
        vertex.append((i*10, i*10+1))
        for i, s in enumerate(sheet[1:]):
            if oseq[i]>0:
                vertex.append((s*10, s*10+1))
            else:
                vertex.append((s*10+1, s*10))
        first = min(vertex)
        if first[0] > first[1]:
            vertex = [(s[1],s[0]) for s in vertex]
        vertices.append(vertex)
    return vertices

def makefatgraph(vertices):
    '''
    Construct vertices, edges, & internal edges from given vertices.
    '''
    vdict = {}
    for vertex in vertices:
        #We want anti-clockwise ordering of half-edges on each vertex
        l = [v[0] for v in vertex] + [v[1] for v in vertex][::-1]
        r = len(vdict)
        for i in range(r, r+len(l)):
            vdict[l[i - r]] = i + 1

    #Now we find edges
    halfedges = sorted([i for vertex in vertices
                        for v in vertex
                        for i in v])
    edges = [(i, j) for i, j in zip(halfedges, halfedges[1:])]
    #edges = [(halfedges[i], halfedges[i+1])
    #         for i in range(1, len(halfedges)-2, 2)]

    #Translate vertices and edges according to vdict
    v = []
    p = 1
    for vertex in vertices:
        v.append(tuple(range(p, 2*len(vertex)+p)))
        p += 2*len(vertex)

    iv = []
    iedges = [e for vertex in vertices for e in vertex]
    for e in iedges:
        iv.append((vdict[e[0]], vdict[e[1]]))

    e = []
    for d in edges:
        edge = (vdict[d[0]], vdict[d[1]])
        if edge[0] > edge[1]:
            edge = edge[::-1]
        if edge not in iv:
            e.append(edge)

    return set(v), set(e), set(iv)

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
    "If s_mat is 2x2, connect the two strands. "
    '''
    array([[ 0,.8,.5,-1, 0],       array([[ 0, 1,-1,-1, 0],
           [ 0, 0,.6,-1,-1],              [ 0, 0, 1,-1,-1],
           [ 0, 0, 0,.4,-1],   |-->       [ 0, 0, 0, 1,-1],
           [ 0, 0, 0, 0,.7],              [ 0, 0, 0, 0, 1],
           [ 0, 0, 0, 0, 0]])             [ 0, 0, 0, 0, 0]])
    '''
    w_mat = np.copy(s_mat)
    m_mat = np.zeros_like(w_mat, dtype=int)
    if len(m_mat) == 2:
        m_mat[0,1] = 1
        return m_mat
    else:
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

def completions(p_mat, o_mat, allow=False):
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
    if allow:
        if not any([isbifurcated(w_mat),
                         hasbarrel(w_mat)]):
            yield np.where(p_mat>0, p_mat, 0), o_mat
    else:
        if not any([isbifurcated(w_mat),
                    hasbarrel(w_mat),
                    hasisolated(w_mat)]):
            yield np.where(p_mat>0, p_mat, 0), o_mat
    for n in range(len(es)//2 + 1):
        if n>0:
            logging.debug('n={}; {}'.format(n-1, len(seen)))
        seen = []
        for indices in combinations_repeat(es, 2, n):
            ps = set(indices)
            if ps in seen:
                continue
            else:
                seen.append(ps)
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
                    hasbarrel(c_mat)]):
                continue
            elif not allow and hasisolated(c_mat):
                continue
            c_mat = np.where(c_mat>0, c_mat, 0)
            #c_mat is valid. Construct all possible omat
            x_mat = np.copy(o_mat)
            for rs in product([False, True], repeat=n):
                for i, v in enumerate(rs):
                    x_mat[indices[i]] = v
                yield np.triu(c_mat, 1), x_mat

def connect(to_pair, new_pair):
    "Connect two sheets with new_pair"
    "to_pair is a set of two sheets"
    '''
    connect({(1,2),(3,4)}, (1,4)) --> (2,1,4,3)
    '''
    if len(to_pair) != 2:
        raise ValueError
    for sheet in to_pair:
        if new_pair[0] == sheet[0]:
            first = sheet[::-1]
        elif new_pair[0] == sheet[-1]:
            first = sheet
        elif new_pair[1] == sheet[0]:
            second = sheet
        elif new_pair[1] == sheet[-1]:
            second = sheet[::-1]
    connected = first + second
    if connected[-1] < connected[0]:
        connected = connected[::-1]
    return connected

def all_pairs(sheets, dist=0, last_pair=None, allow_single=False):
    "Generate all valid pairings of ends, avoiding barrels"
    "Isolated strand is an 1-tuple in sheets."
    '''
    all_pairs({(0,1),(2,3,4)}) -->
    {(0,1), (2,3,4)}, {(0,1,2,3,4)}, {(1,0,2,3,4)}, 
    {(0,1,4,3,2)}, {(2,3,4,0,1)}
    '''
    isolated = set(sheet[0] for sheet in sheets if len(sheet)==1)
    goodsheets = set(sheet for sheet in sheets if len(sheet)>1)
    paired_ends = set(i for sheet in goodsheets
                      for i in (sheet[0], sheet[-1]))

    if allow_single:
        yield sheets
    elif len(isolated) == 0:
        yield sheets
    ends = sorted(list(isolated  | paired_ends))
    for newpair in combinations(ends, 2):
        #newpair must be more than dist apart (-1 in matrix)
        if abs(newpair[0]-newpair[1]) <= dist:
            continue

        #Check for duplicates.
        if last_pair:
            if newpair[0] < last_pair[0]:
                continue
            elif newpair[0]==last_pair[0] and newpair[1]<=last_pair[1]:
                continue

        #Check for barrels
        abandon = False
        for sheet in sheets:
            if newpair[0] in sheet and newpair[1] in sheet:
                abandon = True
                break
        if abandon:
            continue

        #Update sheets
        newsheets = set()
        topair = set()
        for sheet in sheets:
            if set(sheet).isdisjoint(set(newpair)):
                newsheets.add(sheet)
            else:
                topair.add(sheet)
                if len(topair) == 2:
                    newsheet = connect(topair, newpair)
                    newsheets.add(newsheet)
        if len(newsheets) == 1:
            yield newsheets
        else:
            for finalsheets in all_pairs(newsheets, dist,
                                         newpair, allow_single):
                yield finalsheets

def tomat(sheets):
    "Make pairing matrix from sheets"
    n = len(set(i for sheet in sheets for i in sheet))
    mat = np.zeros((n,n), dtype=int)
    for sheet in sheets:
        if len(sheet) < 2:
            continue
        for i, j in zip(sheet, sheet[1:]):
            mat[i,j] = 1
            mat[j,i] = 1
    return mat

def completions2(p_mat, o_mat, allow_single=False):
    "All possible completions of partial pairing matrix, excluding "
    "bifurcations, barrels, and isolated strands."
    "Return a list of tuples consisting score and orientation matrices."
    w_mat = p_mat + p_mat.T
    lst = [idx[1]-idx[0] for idx in zip(*np.where(w_mat<0))]
    try:
        min_dist = max(lst)
    except ValueError:
        min_dist = 0
    zero_mat = np.where(w_mat<0, 0, w_mat)
    sheets = findsheets(zero_mat)
    for goodsheets in all_pairs(sheets, min_dist, None, allow_single):
        c_mat = tomat(goodsheets)
        indices = [idx for idx in zip(*np.nonzero(c_mat-zero_mat))
                   if idx[0] < idx[1]]
        for rs in product([False, True], repeat=len(indices)):
            x_mat = np.copy(o_mat)
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

def reduce_mat(mat):
    "Reduce non-zero entries in mat depending on its size."
    l = len(mat)
    return np.tril(mat, max(l-5, 1))

def make_pmat(pmat, omat):
    "Make single pairing matrix"
    mat = np.zeros_like(pmat)
    for i, j in np.ndindex(pmat.shape):
        if i<j:
            if omat[i,j]:
                mat[i,j] = pmat[i,j]
                mat[j,i] = 0
            else:
                mat[j,i] = pmat[i,j]
                mat[i,j] = 0
    return mat

def repopulate(pmat, original, omat):
    "Re-populate first diagonals in pmat with entries from original "
    "and the rest according to genus-boundary filter."
    with open(MAT_FILE, 'rb') as fh:
        bg_mat = pickle.load(fh)
    sheets = findsheets(pmat)
    vertices = makevertices(sheets, omat)
    v, e, iv = makefatgraph(vertices)
    fg = Fatgraph(v, e)
    g = fg.genus
    b = len(fg.boundaries)
    k = max(len(sheet) for sheet in sheets)
    try:
        score = bg_mat[g,b,k] / np.sum(bg_mat[:,:,k])
    except IndexError:
        score = 0
    l = len(pmat)
    n = min(5, max(l-5, 1))
    rmat = np.triu(pmat*score, n+1) + np.tril(original, n)
    return np.where(rmat>0, rmat, 0)

def main(partial, allow, save, dbg):
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

    part_mat = makepartialmatrix(s_mat)
    part_mat = reduce_mat(part_mat)
    logger.debug('Partial matrix:')
    logger.debug(part_mat)
    sum_mat = np.zeros_like(s_mat)
    cnt = 0

    for pmat, omat in completions2(part_mat, o_mat, allow):
        pmat = repopulate(pmat, s_mat, omat)
        mat = make_pmat(pmat, omat)
        sum_mat += mat
        cnt += 1
    logger.debug('{} completions found.'.format(cnt))
    mat = sum_mat / cnt
    if save:
        sd = pathlib.Path(save)
        of = sd / '{}.pkl'.format(pid)
        with open(of, 'wb') as fh:
            pickle.dump(mat, fh)
    else:
        print(mat)

    end = datetime.now()
    d = end - start
    logger.info('{}: Finished {} in {} seconds'.format(
        end, pid, d.seconds))

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Construct protein metastructure from partial '
        'motif by going through all possible completions of partial '
        'matrix and applying genus filter.')
    parser.add_argument('partial',
                        help='Partial motif pkl file.')
    parser.add_argument('-i', '--allow-isolated', action='store_true',
                        help='Allow isolated strand in the resulting '
                        'motifs.')
    parser.add_argument('-s', '--save',
                        help='Save output to this directory.')
    parser.add_argument('-d', '--debug', action='store_true',
                        help='Print debug messages.')
    args = parser.parse_args()
    main(args.partial, args.allow_isolated, args.save, args.debug)
