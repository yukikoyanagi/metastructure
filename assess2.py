#!/usr/bin/env python3

MAX_SIZE = 20
MAT_FILE = 'gbs_arr_20.pkl'
#MAT_FILE = 'gbs_nstr_20.pkl'
#MAT_FILE = 'genbound_arr.pkl'
BG_CUTOFF = 0.01

import argparse, pickle, pathlib, logging
import multiprocessing as mp
from itertools import accumulate
from operator import mul
import numpy as np

from fatgraph.fatgraph import Fatgraph

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
    --> [[0,2,1],[3]]
    '''
    wmat = np.triu(p_mat,1) + np.triu(p_mat,1).T
    edges = [i for i in range(len(wmat))
             if np.count_nonzero(wmat[i])<2]
    sheets = []
    while edges:
        e = edges.pop(0)
        sheet = findsheetat(e, wmat)
        if sheet[0] > sheet[-1]:
            sheet = sheet[::-1]
        sheets.append(sheet)
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

def compare(tmat, cmat, ori=False):
    if not ori:
        tmat = np.triu(tmat + tmat.T)
        cmat = np.triu(cmat + cmat.T)
    TP = np.count_nonzero(tmat*cmat)
    FP = np.count_nonzero(cmat-tmat>0)
    FN = np.count_nonzero(tmat-cmat>0)
    return TP, FP, FN

def mot2mat(mot, l):
    "Make complete pairing matrix from motif"
    pmat = np.zeros((l,l), dtype=int)
    for link in mot:
        i, j = link[0]
        if i > j:
            i, j = j, i
        if not link[1]:
            i, j = j, i
        pmat[i,j] = 1
    return pmat

def makepomat(mat):
    "Convert single pairing matrix to pairing and orientation matrices"
    pmat = np.triu(mat + mat.T)
    omat = np.triu(mat)
    return pmat, omat

def applyfilter(mat, bg_mat):
    logger = logging.getLogger()
    pmat, omat = makepomat(mat)
    sheets = findsheets(pmat)
    vertices = makevertices(sheets, omat)
    v, e, iv = makefatgraph(vertices)
    fg = Fatgraph(v, e)
    g = fg.genus
    b = len(fg.boundaries)
    k = max(len(sheet) for sheet in sheets) #max. sheet size
    #k = sum(len(sheet) for sheet in sheets) #sheet count
    if bg_mat.ndim == 2:
        try:
            score = bg_mat[g,b] / np.sum(bg_mat)
        except IndexError:
            score = 0
    elif bg_mat.ndim == 3:
        try:
            score = bg_mat[g,b,k] / np.sum(bg_mat[:,:,k])
        except IndexError:
            score = 0
    return bool(score >= BG_CUTOFF) #return bool, not numpy.bool_

def compute(tf, pf, size, orientation=False):
    "Compute filtering results for a single candidate file." 
    pid = tf.stem
    with open(MAT_FILE, 'rb') as fh:
        bg_mat = pickle.load(fh)

    output = []

    with open(pf, 'rb') as fh:
        cmots = pickle.load(fh)

    with open(tf, 'rb') as fh:
        themot = pickle.load(fh)
    tmot = set()
    for link in themot:
        pair, ori = link
        if pair[0] < pair[1]:
            tmot.add(((pair[0], pair[1]), not ori[0]==ori[1]))
        else:
            tmot.add(((pair[1], pair[0]), not ori[0]==ori[1]))
    try:
        tmat = mot2mat(tmot, size)
    except IndexError as e:
        print(tf.name, tmot)
        raise e
    
    for cmot in cmots:
        try:
            cmat = mot2mat(cmot, size)
        except IndexError as e:
            print(tf.name, cmot)
            raise e
        accepted = applyfilter(cmat, bg_mat)
        TP, FP, FN = compare(tmat, cmat, orientation)
        TPR = TP / (TP + FN)
        try:
            PPV = TP / (TP + FP)
        except ZeroDivisionError:
            PPV = 0.0
        output.append((pid, size, accepted, TPR, PPV))
    return output

def main(tdir, pdir, nstr, v, orientation, cpus, msize, output):
    logging.basicConfig(level=logging.INFO)

    global BG_CUTOFF
    BG_CUTOFF = v
    global MAX_SIZE
    MAX_SIZE = msize

    pdir = pathlib.Path(pdir)
    tdir = pathlib.Path(tdir)

    sizes = dict()
    with open(nstr) as fh:
        for line in fh:
            sizes[line.split()[0]] = int(line.split()[1])

    inputs = []
    for mf in pdir.glob('*.pkl'):
        if len(mf.name) > 9:
            fn = '{}.pkl'.format(mf.name[:5])
        else:
            fn = mf.name
        tf = tdir / '{}'.format(fn)
        size = sizes[tf.stem]
        if size > MAX_SIZE:
            continue

        inputs.append((tf, mf, size))
    pool = mp.Pool(processes=cpus)
    res_objects = [pool.apply_async(
        compute, args=(i[0], i[1], i[2], orientation)
    ) for i in inputs]
    pool.close()
    pool.join()
    results = [item for r in res_objects for item in r.get()]

    if output:
        with open(output, 'wb') as fh:
            pickle.dump(results, fh, pickle.HIGHEST_PROTOCOL)
    else:
        print(results)



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Assess candidate  '
                                     'metastructures using genus-'
                                     'boundary filter.')
    parser.add_argument('true_motif_dir',
                        help='Directory containing pickle files '
                        'of true motifs.')
    parser.add_argument('pred_motif_dir',
                        help='Directory containing pickle files '
                        'of candidate motifs.')
    parser.add_argument('n_strands',
                        help='File containing number of strands per '
                        'protein.')
    parser.add_argument('-v', '--cutoff', type=float, default=BG_CUTOFF,
                        help='Cutoff value for heatmap filter')
    parser.add_argument('-o', '--orientation', action='store_true',
                        help='Consider parallel/anti-parallel when '
                        'computing accuracy.')
    parser.add_argument('-c', '--cpus', type=int,
                        help='Number of cores to use in computation')
    parser.add_argument('-m', '--maxsize', type=int, default=MAX_SIZE,
                        help='Maximum size of protein to consider')
    parser.add_argument('-s', '--save',
                        help='Save results in this file.')
    args = parser.parse_args()
    main(args.true_motif_dir, args.pred_motif_dir, args.n_strands,
         args.cutoff, args.orientation, args.cpus, args.maxsize,
         args.save)

