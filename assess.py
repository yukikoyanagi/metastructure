#!/usr/bin/env python3

MAX_SIZE = 11
MAT_FILE = 'gbs_arr.pkl'
#MAT_FILE = 'genbound_arr.pkl'
BG_CUTOFF = 0.01

import argparse, pickle, pathlib
from itertools import combinations, accumulate
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


def compare(tmot, pmot, ori):
    tpairs = {t[0] for t in tmot}
    ppairs = {p[0] for p in pmot}
    tlinks = {t for t in tmot}
    plinks = {p for p in pmot}
    strands = range(max(s for p in tpairs for s in p) + 1)
    apairs = set(combinations(strands, 2))
    if ori:
        TP = len(tlinks & plinks)
        FP = len(plinks - tlinks)
    else:
        TP = len(tpairs & ppairs)
        FP = len(ppairs - tpairs)
    TN = len((apairs - tpairs) & (apairs - ppairs))
    FN = len(tpairs - ppairs)

    c_pairs = set(p for p in pmot if p[0] in tpairs)
    c_links = set(p for p in pmot if p in tmot and p[0] in tpairs)

    return TP, FP, TN, FN, len(c_pairs), len(c_links)

def write_assess(tmots, cmots, ori):
    a = {}
    for k in tmots:
        mots = cmots[k]
        tmot = tmots[k]
        acus = []
        for pmot, acc in mots:
            TP, FP, TN, FN, CP, CL = compare(tmot, pmot, ori)
            acu = (TP + TN) / (TP + FP + TN + FN)
            acus.append((acu, acc))
        a[k] = acus

    l = []
    hi_acc = 0
    hi_rej = 0
    for k in a:
        try:
            hi_acc = max(p for p,t in a[k] if t)
        except ValueError:
            hi_acc = 0
        try:
            hi_rej = max(p for p,t in a[k] if not t)
        except ValueError:
            hi_rej = 0
        l.append(hi_acc >= hi_rej)
    print(
        'Highest accuray accepted: {}/{} ({:.2%})'.format(
            sum(l), len(l), sum(l)/len(l))
    )
    
    n = sum(len(v) for v in a.values())
    m = len([w for v in a.values() for w, t in v if w>=0.8])
    print(
        '# of candidates at 80% accuracy level: {:>5}/{:>5} ({:>6.2%})'.format(
            m, n, m/n)
    )
    h = sum(t for v in a.values() for w, t in v if w>=0.8)
    print(
        '# of acceptance at 80% accuracy level: {:>5}/{:>5} ({:>6.2%})'.format(
            h, m, h/m)
    )
    m = len([w for v in a.values() for w, t in v if w>=0.9])
    print(
        '# of candidates at 90% accuracy level: {:>5}/{:>5} ({:>6.2%})'.format(
            m, n, m/n)
    )
    h = sum(t for v in a.values() for w, t in v if w>=0.9)
    print(
        '# of acceptance at 90% accuracy level: {:>5}/{:>5} ({:>6.2%})'.format(
            h, m, h/m)
    )
    m = len([w for v in a.values() for w, t in v if w>=0.95])
    print(
        '# of candidates at 95% accuracy level: {:>5}/{:>5} ({:>6.2%})'.format(
            m, n, m/n)
    )
    h = sum(t for v in a.values() for w, t in v if w>=0.95)
    print(
        '# of acceptance at 95% accuracy level: {:>5}/{:>5} ({:>6.2%})'.format(
            h, m, h/m)
    )
    m = len([w for v in a.values() for w, t in v if w>=0.99])
    print(
        '# of candidates at 99% accuracy level: {:>5}/{:>5} ({:>6.2%})'.format(
            m, n, m/n)
    )
    h = sum(t for v in a.values() for w, t in v if w>=0.99)
    print(
        '# of acceptance at 99% accuracy level: {:>5}/{:>5} ({:>6.2%})'.format(
            h, m, h/m)
    )
    print('\n')

def write_accepted(tmots, cmots):
    with open(MAT_FILE, 'rb') as fh:
        bg_mat = pickle.load(fh)
    dim = bg_mat.ndim
    n = MAX_SIZE + 1
    accepted = [0 for i in range(n)]
    rejected = [0 for i in range(n)]

    for k in tmots:
        mots = cmots[k]
        tmot = tmots[k]
        l = len(set(i for p, o in tmot for i in p))
        for mot in mots:
            if mot[-1]:
                accepted[l] += 1
            else:
                rejected[l] += 1

    print('==========================================================')
    print('\n')
    print('Number of candidates accepted at v={} ({}D heatmap filter)'.format(BG_CUTOFF, dim))
    print('# of strands:' + ''.join(['{:>8}'.format(i) for i in range(n)]))
    print('Accepted    :' + ''.join(['{:>8}'.format(accepted[i]) for i in range(n)]))
    print('Rejected    :' + ''.join(['{:>8}'.format(rejected[i]) for i in range(n)]))
    pcts = [0 for i in range(n)]
    for i in range(n):
        try:
            pcts[i] = accepted[i] / (accepted[i] + rejected[i])
        except ZeroDivisionError:
            pcts[i] = 0
    print('Accepted %  :' + ''.join(['{:>8.1%}'.format(pcts[i]) for i in range(n)]))
    print('\nTotal accepted: {:.2%}'.format(sum(accepted) / (sum(accepted) + sum(rejected)) ))

def mot2mat(mot):
    "Make pairing and orientation matrices from motif"
    l = len(set(i for p, o in mot for i in p))
    p_mat = np.zeros((l,l), dtype=int)
    o_mat = np.zeros((l,l), dtype=bool)
    for link in mot:
        i, j = link[0]
        p_mat[i,j] = 1
        o_mat[i,j] = link[1]
    return p_mat, o_mat

def apply_filter(mots, v):
    global BG_CUTOFF
    BG_CUTOFF = v
    with open(MAT_FILE, 'rb') as fh:
        bg_mat = pickle.load(fh)
    cmots = []
    for mot in mots:
        pmat, omat = mot2mat(mot)
        sheets = findsheets(pmat)
        vertices = makevertices(sheets, omat)
        v, e, iv = makefatgraph(vertices)
        fg = Fatgraph(v, e)
        g = fg.genus
        b = len(fg.boundaries)
        k = max(len(sheet) for sheet in sheets)
        if bg_mat.ndim == 2:
            try:
                score = bg_mat[g,b] / np.sum(bg_mat)
            except IndexError:
                score = 0
        elif bg_mat.ndim == 3:
            try:
                score = bg_mat[g,b,k] / np.sum(bg_mat)
            except IndexError:
                score = 0
        cmots.append((mot, score >= BG_CUTOFF))
    return cmots

def main(tdir, pdir, v, ori):

    pdir = pathlib.Path(pdir)
    tdir = pathlib.Path(tdir)

    cmots = {}
    tmots = {}
    for mf in pdir.glob('*.pkl'):
        with open(mf, 'rb') as fh:
            mots = pickle.load(fh)
        cmots[mf.stem] = mots

        with open(tdir / '{}'.format(mf.name), 'rb') as fh:
            themot = pickle.load(fh)

        tmot = set()
        for link in themot:
            pair, ori = link
            if pair[0] < pair[1]:
                tmot.add(((pair[0], pair[1]), not ori[0]==ori[1]))
            else:
                tmot.add(((pair[1], pair[0]), not ori[0]==ori[1]))
        tmots[mf.stem] = tmot

    for k in cmots:
        cmots[k] = apply_filter(cmots[k], v)
    write_accepted(tmots, cmots)
    write_assess(tmots, cmots, ori)



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Assessment of '
                                     'genus-boundary filter for '
                                     'selecting good candidate '
                                     'structures.')
    parser.add_argument('true_motif_dir',
                        help='Directory containing pickle files '
                        'of true motifs.')
    parser.add_argument('pred_motif_dir',
                        help='Directory containing pickle files '
                        'of candidate motifs.')
    parser.add_argument('-v', '--cutoff', type=float, default=0.01,
                        help='Cutoff value for heatmap filter')
    parser.add_argument('-o', '--orientation', action='store_true',
                        help='Consider parallel/anti-parallel when '
                        'computing accuracy.')
    args = parser.parse_args()
    main(args.true_motif_dir, args.pred_motif_dir,
         args.cutoff, args.orientation)
