#!/usr/bin/env python3

MAX_SIZE = 18
MAT_FILE = 'gbs_arr.pkl'
#MAT_FILE = 'genbound_arr.pkl'
BG_CUTOFF = 0.01

import argparse, pickle, pathlib, logging
import multiprocessing as mp
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


def compare(tmot, pmot, ori=False):
    tpairs = {t[0] for t in tmot}
    ppairs = {p[0] for p in pmot}
    tlinks = {t for t in tmot}
    plinks = {p for p in pmot}
    strands = range(max(s for p in tpairs for s in p) + 1)
    apairs = set(combinations(strands, 2))
    if ori:
        l = len(strands)
        tmat = mot2mat2(tmot, l)
        pmat = mot2mat2(pmot, l)
        TP = np.count_nonzero(tmat*pmat)
        FP = np.count_nonzero(pmat-tmat>0)
        TN = np.count_nonzero(pmat+tmat==0) - l
        FN = np.count_nonzero(tmat-pmat>0)
    else:
        TP = len(tpairs & ppairs)
        FP = len(ppairs - tpairs)
        TN = len((apairs - tpairs) & (apairs - ppairs))
        FN = len(tpairs - ppairs)

    c_pairs = set(p for p in pmot if p[0] in tpairs)
    c_links = set(p for p in pmot if p in tmot and p[0] in tpairs)

    return TP, FP, TN, FN, len(c_pairs), len(c_links)

def write_accepted(results, by_protein, count_proteins):
    #vs = [.8, .85, .9, .95, 1]
    vs = [.5, .6, .7, .8, .9, 1]
    byvar_dic = {3: 'accuracy', 4: 'sensitivity', 5: 'specificity', 6: 'precision'}
    byvar = 4 #Index of var to measure the result. 3=accuracy
    with open(MAT_FILE, 'rb') as fh:
        bg_mat = pickle.load(fh)
    dim = bg_mat.ndim
    n = MAX_SIZE + 1
    nall = len(results)
    accepted = [0 for i in range(n)]
    rejected = [0 for i in range(n)]
    nc = [0 for _ in vs] #number of candidates
    ac = [0 for _ in vs] #accepted candidates
    if by_protein:
        np = [0 for _ in vs] #number of proteins
        ap = [0 for _ in vs] #accepted proteins

    l = []
    hi_acc = 0
    hi_rej = 0
    n_cand = [[] for i in range(n)]
    cand_count = 0
    pid = ''

    for r in results:
        if pid and r[0] != pid:
            #Process one pid
            l.append(hi_acc >= hi_rej)
            hi_acc = 0
            hi_rej = 0
            n_cand[last_size].append(cand_count)
            cand_count = 0
            if by_protein: #add np to nc
                nc = [x+y for x,y in zip(nc, np)]
                ac = [x+y for x,y in zip(ac, ap)]
                np = [0 for _ in vs]
                ap = [0 for _ in vs]
        pid = r[0]
        last_size = r[1]

        #Highest accuracy among accepted/rejected
        if r[2]:
            accepted[r[1]] += 1
            if r[byvar] > hi_acc:
                hi_acc = r[byvar]
        else:
            rejected[r[1]] += 1
            if r[byvar] > hi_rej:
                hi_rej = r[byvar]

        #Acceptance by minimum accuracy
        for i, v in enumerate(vs):
            if r[byvar] >= v:
                if by_protein:
                    np[i] = 1
                else:
                    nc[i] += 1
                if r[2]:
                    if by_protein:
                        ap[i] = 1
                    else:
                        ac[i] += 1

        #Count candidates
        cand_count += 1

    #Still need to process last pid
    l.append(hi_acc >= hi_rej)
    n_cand[r[1]].append(cand_count)
    if by_protein:
        nc = [x+y for x,y in zip(nc, np)]
        ac = [x+y for x,y in zip(ac, ap)]

    print('==========================================================')
    print('\n')
    if count_proteins:
        
        print('Number of proteins')
        print('# Strands  :' + ''.join(['{:>8}'.format(i) for i in range(n)]))
        print('# Proteins :' + ''.join(['{:>8}'.format(len(c)) for c in n_cand]))
        print('#cands max.:' + ''.join(['{:>8}'.format(max(c, default=0)) for c in n_cand]))
        print('#cands min.:' + ''.join(['{:>8}'.format(min(c, default=0)) for c in n_cand]))
        avgs = []
        for c in n_cand:
            if len(c):
                avgs.append(sum(c)/len(c))
            else:
                avgs.append(0)
        print('Avg. cands :' + ''.join(['{:>8.1f}'.format(a) for a in avgs]))
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
    by_var = byvar_dic[byvar]
    print(
        'Highest {} accepted: {}/{} ({:.2%})'.format(
            by_var, sum(l), len(l), sum(l)/len(l))
    )
    if by_protein:
        c = ' proteins '
        nall = len(l)
    else:
        c = 'candidates'

    for i, v in enumerate(vs):
        print(
            '# of {} at {:>4.0%} {} level: {:>8}/{:>8} ({:>6.2%})'.format(
                c, v, by_var, nc[i], nall, nc[i]/nall
            )
        )
        print(
            '# of acceptance at {:>4.0%} {} level: {:>8}/{:>8} ({:>6.2%})'.format(
                v, by_var, ac[i], nc[i], ac[i]/nc[i]
            )
        )
    print('\n')

def mot2mat(mot, l):
    "Make pairing and orientation matrices from motif"
    p_mat = np.zeros((l,l), dtype=int)
    o_mat = np.zeros((l,l), dtype=bool)
    for link in mot:
        i, j = link[0]
        p_mat[i,j] = 1
        o_mat[i,j] = link[1]
    return p_mat, o_mat

def mot2mat2(mot, l):
    "Make complete pairing matrix from motif"
    pmat = np.zeros((l,l), dtype=int)
    for link in mot:
        i, j = link[0]
        if not link[1]:
            i, j = j, i
        pmat[i,j] = 1
    return pmat

def apply_filter(mots, v, size):

    global BG_CUTOFF
    BG_CUTOFF = v
    with open(MAT_FILE, 'rb') as fh:
        bg_mat = pickle.load(fh)
    cmots = []
    for mot in mots:
        pmat, omat = mot2mat(mot, size)
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
                score = bg_mat[g,b,k] / np.sum(bg_mat[:,:,k])
            except IndexError:
                score = 0
        cmots.append((mot, score >= BG_CUTOFF))
    return cmots

def applyfilter(mot, bg_mat, size):
    logger = logging.getLogger()
    try:
        pmat, omat = mot2mat(mot, size)
    except IndexError as e:
        logger.info('{} with size {}'.format(mot, size))
        raise e
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
            score = bg_mat[g,b,k] / np.sum(bg_mat[:,:,k])
        except IndexError:
            score = 0
    return score >= BG_CUTOFF

def computemany(data, o=False):
    logger = logging.getLogger()
    output = []
    for pid, tmot, cmot in data:
        try:
            TP, FP, TN, FN, CP, CL = compare(tmot, cmot, o)
        except IndexError as e:
            print('IndexError in {}'.format(pid))
            raise e
        accuracy = (TP + TN) / (TP + TN + FP + FN)
        sensitivity = TP / (TP + FN)
        try:
            specificity = TN / (TN + FP)
        except ZeroDivisionError as e: #this could happen in 2-strand
            specificity = 1.0
        try:
            precision = TP / (TP + FP)
        except ZeroDivisionError as e:
            precision = 0.0
        size = len(set(e for p, _ in tmot for e in p))
        with open(MAT_FILE, 'rb') as fh:
            bg_mat = pickle.load(fh)
        accepted = applyfilter(cmot, bg_mat, size)
        output.append((pid, size, accepted, accuracy,
                       sensitivity, specificity, precision))
    return output

def main(tdir, pdir, v, orientation, cpus, msize, cnt, byprot, raw):
    logging.basicConfig(level=logging.INFO)

    global BG_CUTOFF
    BG_CUTOFF = v
    global MAX_SIZE
    MAX_SIZE = msize

    pdir = pathlib.Path(pdir)
    tdir = pathlib.Path(tdir)

    cmots = {}
    tmots = {}
    for mf in pdir.glob('*.pkl'):
        with open(tdir / '{}'.format(mf.name), 'rb') as fh:
            themot = pickle.load(fh)
        if max(e for p, _ in themot for e in p) > MAX_SIZE - 1:
            continue

        tmot = set()
        for link in themot:
            pair, ori = link
            if pair[0] < pair[1]:
                tmot.add(((pair[0], pair[1]), not ori[0]==ori[1]))
            else:
                tmot.add(((pair[1], pair[0]), not ori[0]==ori[1]))
        tmots[mf.stem] = tmot

        with open(mf, 'rb') as fh:
            mots = pickle.load(fh)
        cmots[mf.stem] = mots

    i = 0
    data = []
    chunk = []
    for k in cmots:
        for mot in cmots[k]:
            chunk.append((k, tmots[k], mot))
            i += 1
            if i == 5000:
                data.append(chunk)
                i = 0
                chunk = []
    else:
        data.append(chunk)
    
    pool = mp.Pool(processes=cpus)
    res_objects = [pool.apply_async(
        computemany, args=(chunk, orientation)
    ) for chunk in data]
    pool.close()
    pool.join()
    results = [item for r in res_objects for item in r.get()]

    if raw:
        with open(raw, 'wb') as fh:
            pickle.dump(results, fh)
    else:
        write_accepted(results, byprot, cnt)



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
    parser.add_argument('-c', '--cpus', type=int,
                        help='Number of cores to use in computation')
    parser.add_argument('-m', '--maxsize', type=int,
                        help='Maximum size of protein to consider')
    parser.add_argument('-p', '--count_proteins', action='store_true',
                        help='Display count of proteins per size')
    parser.add_argument('-b', '--by_protein', action='store_true',
                        help='Produce acceptance results by protein')
    parser.add_argument('-r', '--raw_results',
                        help='Save raw results in this file and exit')
    args = parser.parse_args()
    main(args.true_motif_dir, args.pred_motif_dir,
         args.cutoff, args.orientation, args.cpus, args.maxsize,
         args.count_proteins, args.by_protein, args.raw_results)
