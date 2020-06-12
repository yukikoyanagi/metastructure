#!/usr/bin/env python3

SUM_DIR = '/home/au447708/QGM/metastr/data/summaries2020'
#SUM_DIR = '/home/au447708/QGM/metastr/data/summaries_cov2'
HQ60_DIR = '/home/au447708/QGM/cdp/data/HQ60_20'
#HQ60_DIR = '/home/au447708/QGM/metastr/data/HQ60_cov2'

#column numbers in summary* file
RESNUM_COL = 1
DSSP_COL = 3
DON_COL = 10
ACC_COL = 11


import argparse, re, pickle
import numpy
from operator import attrgetter, itemgetter, mul
from itertools import accumulate, combinations

def getseq(name):
    '''
    Get sequence of 3-class DSSP from protein name (id & chain).
    DSSP chars are ordered by DSSP residue numbering. 
    In addition to H, S, C, the resulting sequence may contain 
    ! for chain break, and/or previous chains (chains B, C, etc. 
    may start at residue number other than 1).
    '''
    seq = []
    with open('{}/summary{}.txt'.format(SUM_DIR, name.upper())) as fh:
        for line in fh:
            cols = line.split()
            n = int(cols[RESNUM_COL])
            c = cols[DSSP_COL]
            if n-1 == len(seq):
                seq.append(c)
            elif n-1 > len(seq):
                seq.extend(['!']*(n-1-len(seq)))
                seq.append(c)
            else:
                seq[n-1] = c
    # Single ! sandwiched by two S's is also an S
    for i in range(len(seq)):
        if seq[i] == '!' and seq[i-1] == 'S' and seq[i+1] == 'S':
            seq[i] = 'S'
        # ! franked by an S is also an S
        try:
            if seq[i] == '!' and (
                    (seq[i-1] == 'S' and seq[i-2] != 'S') or
                    (seq[i+1] == 'S' and seq[i+2] != 'S')):
                seq[i] = 'S'
        except IndexError:
            continue

    return ''.join(seq)

def getbonds(name):
    '''
    Return a list of bonds [(from, to),], with each bond expressed
    as a 2-tuple of DSSP residue numbers. 
    '''
    bonds = []
    p = re.compile(r'__[US]{2}')
    with open('{}/{}00.txt'.format(HQ60_DIR, name.upper())) as fh:
        for line in fh:
            if p.search(line):
                cols = line.split()
                don = int(cols[DON_COL])
                acc = int(cols[ACC_COL])
                bonds.append((don, acc))
    return bonds

def getidx(mat, row):
    '''
    Return index of row in matrix. If there are more than one 
    matching rows, return the index of the first match.
    '''
    return numpy.where(numpy.all(mat==row, axis=1))[0][0].item()

def extract(strands, mat, ridx, bifurcation=False, barrel=False):
    '''
    Extract a vertex from strands. The first strand is given by the
    corresponding row in mat.
    '''
    vertex = []
    idx = ridx
    #print('idx is {}'.format(idx))
    vertex.append(strands[idx])
    try:
        nidx = numpy.flatnonzero(mat[idx]).item()
    except ValueError:
        nidx = numpy.flatnonzero(mat.T[idx]).item()
    #print('nidx is {}'.format(nidx))
    vertex.append(strands[nidx])
    nrow = mat[nidx]
    crow = numpy.copy(nrow)
    crow[idx] = 0
    ccol = numpy.copy(mat.T[nidx])
    ccol[idx] = 0
    while (numpy.count_nonzero(crow) or
           numpy.count_nonzero(ccol)):
        if (numpy.count_nonzero(crow) > 1 or
            numpy.count_nonzero(ccol) >1):
            #we have a bifurcation
            if bifurcation:
                raise ValueError('Bifurcation detected.')
            exit()
        try:
            idx = numpy.flatnonzero(crow).item()
        except ValueError:
            idx = numpy.flatnonzero(ccol).item()
        if strands[idx] in vertex:
            # we have a barrel
            if barrel:
                raise ValueError('Beta-barrel detected.')
            exit()
        vertex.append(strands[idx])
        crow = numpy.copy(mat[idx])
        crow[nidx] = 0
        ccol = numpy.copy(mat.T[idx])
        ccol[nidx] = 0
        nidx = idx
        #print('nidx is {}'.format(nidx))
    return vertex

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
            vdict[l[i - r]] = i+1

    #Now we find edges
    halfedges = sorted([i for vertex in vertices
                        for v in vertex
                        for i in v])
    edges = [(halfedges[i], halfedges[i+1])
             for i in range(1, len(halfedges)-2, 2)]

    #Translate vertices and edges according to vdict
    v = []
    p = 0
    for vertex in vertices:
        v.append(tuple(range(p+1, 2*len(vertex)+p+1)))
        p += 2*len(vertex)

    iv = []
    iedges = [e for vertex in vertices for e in vertex]
    for e in iedges:
        iv.append((vdict[e[0]], vdict[e[1]]))

    e = []
    for d in edges:
        e.append((vdict[d[0]], vdict[d[1]]))

    return v, e, iv

def isparallel(f, s, bonds, t=1):
    '''
    True if the two strands f and s are parallel
    '''
    ftos = [b for b in bonds
            if b[0] in range(min(f), max(f)+1)
            and b[1] in range(min(s), max(s)+1)]
    if len(ftos) > 1:
        ftos = sorted(ftos, key=itemgetter(0))
        return ftos[0][1] < ftos[-1][1]
    stof = [b for b in bonds
            if b[0] in range(min(s), max(s)+1)
            and b[1] in range(min(f), max(f)+1)]
    if len(stof) > 1:
        stof = sorted(stof, key=itemgetter(0))
        return stof[0][1] < stof[-1][1]
    if ftos and stof: # neither is empty
        return not ftos[0] == (stof[0][1], stof[0][0])
    # We only have one bond to determin orientation.
    # Try extending strands by up to three residues on either side.
    # If it is still not possible to determine orientation,
    # assume parallel. 
    # We assume there would be two bonds if it is anti-parallel
    if t < 3:
        newf = (f[0]-1, f[1]+1)
        news = (s[0]-1, s[1]+1)
        return isparallel(newf, news, bonds, t+1)
    if t == 3:
        return True
    

def getplacement(start, end):
    if start[0] < start[-1]:
        first = 'r'
    else:
        first = 'l'
    if end[0] < end[-1]:
        second = 'l'
    else:
        second = 'r'
    return first + second

def getmotif(vertices):
    '''
    Return motif for the metastructure graph. A motif is a collection
    of chords which represent parallel h-bonds between strands in a 
    beta-sheet. Each chord is represented by a 2-tuple:
    1st element is a 2-tuple of int's showing the start and end 
    strands of a chord (always start < end). 2nd element is one of
    four strings 'll', 'lr', 'rl', 'rr', representing the attachment
    of the chord. 'll' means the chord is attached on the left side
    of start-strand and end-strand w.r.t. the orientation of the 
    backbone, and so on.
    '''
    rmotif = []
    for vertex in vertices:
        for start, end in zip(vertex, vertex[1:]):
            sign = getplacement(start, end)
            rmotif.append(((start, end), sign))

    sdict = {}
    strands = sorted([s for v in vertices for s in v])
    for i in range(len(strands)):
        sdict[strands[i]] = i

    motif = []
    for rchord in rmotif:
        s, e = rchord[0]
        chord = ((sdict[s], sdict[e]), rchord[1])
        motif.append(chord)
    return set(motif)

def getmotif2(vertices):
    "Return motif in a form {(pair, parallel)}"
    # getmotif2([[(10,12), (22,20), (17,19)]]) --> {((1,2), False), ((0,2), False)}
    mot = getmotif(vertices)
    motif = set()
    for p in mot:
        pair = p[0]
        if pair[0] > pair[1]:
            pair = (pair[1], pair[0])
        o = p[1]
        if o[0] == o[1]:
            parallel = False
        else:
            parallel = True
        motif.add((pair, parallel))
    return motif

def main(name, bif, bar, mot, std, save, debug):
    try:
        seq = getseq(name)
    except FileNotFoundError:
        exit()
    try:
        bonds = getbonds(name)
    except FileNotFoundError:
        exit()

    if debug:
        print(name, seq)
    strands = []
    for m in re.finditer(r'S+', seq):
        strands.append((m.start()+1, m.end()))

    if std:
        print('{}:{}'.format(name, len(strands)))
        exit()

    #Construct nxn matrix, where n=len(strands). (i,j)-entry shows
    #the number of bonds from i'th strand to j'th strand.
    #Also discard the bonds which are not between strands.
    mat = numpy.zeros((len(strands), len(strands)), dtype=int)
    lstr = [list(range(s[0], s[1]+1)) for s in strands]
    newbonds = []
    for bond in bonds:
        don = acc = None
        for s in lstr:
            if bond[0] in s:
                don = lstr.index(s)
            elif bond[1] in s:
                acc = lstr.index(s)
        if don is not None and acc is not None:
            mat[don, acc] += 1
            newbonds.append(bond)
    sbonds = newbonds

    if debug:
        print(strands)
        print(bonds)
        print(sbonds)
        print(mat)

    #Construct vertices
    # isol is used to keep track of the # of isolated strands, which
    # we assume to be the results of interchain bonds; the isolated
    # strand is in fact bonded to another chain.
    vertices = []
    isol = 0
    for idx in range(len(mat)):
        row = mat[idx]
        col = mat.T[idx]
        if (numpy.count_nonzero(row) + numpy.count_nonzero(col) == 0):
            isol += 1
        if ((numpy.count_nonzero(row) + numpy.count_nonzero(col) == 1) or
            # row and col have only one nonzero entry together
            (numpy.count_nonzero(row) ==1 and
             numpy.count_nonzero(col) ==1 and
             numpy.flatnonzero(row)[0] == numpy.flatnonzero(col)[0])
            # row and col have one nonzero entry each, but they are
            # for the same strand
        ):
            if strands[idx] in [s for v in vertices for s in v]:
                continue
            try:
                vertices.append(extract(strands, mat, idx, bif, bar))
            except ValueError as e:
                print('{}: {}'.format(name, e))
                exit()

    # If we don't have all strands in vertices here, we have a barrel
    if len([s for v in vertices for s in v]) < len(strands) - isol:
        if bar:
            print('{}: Barrel detected.'.format(name))
        exit()
    # Make sure the vertices do not contain duplicates
    nvs = []
    for v in vertices:
        if all([(not set(v) < set(w)) for w in vertices]):
            nvs.append(v)
    vertices = nvs

    if debug:
        print(vertices)

    if bif or bar:
        # We only indicate bifurcation or barrel
        exit()

    # Now there is a choice of two orderings of strands in each vertex.
    # We pick the one where the top strand comes before the bottom
    # strands in backbone.
    nvs = []
    for v in vertices:
        if v[0] < v[-1]:
            nvs.append(v)
        else:
            v.reverse()
            nvs.append(v)
    vertices = nvs

    # Determine parallel/anti-parallel
    nvs = []
    for v in vertices:
        parallel_seq = []
        for first, second in zip(v, v[1:]):
            parallel_seq.append(isparallel(first, second, bonds))
        pseq = [1 if i else -1 for i in parallel_seq]
        pseq = list(accumulate(pseq, mul))
        newv = [v[0],]
        for i in range(1,len(v)):
            if pseq[i-1] == 1:
                newv.append(v[i])
            else:
                newv.append((v[i][1], v[i][0]))
        nvs.append(newv)
    vertices = nvs

    if debug:
        print(vertices)

    if save and mot:
        with open('{}/{}.pkl'.format(save.rstrip('/'), name), 'wb') as fh:
            pickle.dump(getmotif(vertices), fh)
        exit()
    if mot:
        print(getmotif2(vertices))
        exit()

    vertices, edges, in_edges = makefatgraph(vertices)

    if save:
        with open('{}/{}.pkl'.format(save.rstrip('/'), name), 'wb') as fh:
            pickle.dump((vertices, edges, in_edges, seq), fh)
    else:
        print('{}:vertices={}'.format(name, vertices))
        print('{}:edges={}'.format(name, edges))
        print('{}:in_edges={}'.format(name, in_edges))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('name',
                        help="Protein's name, including chain id.")
    parser.add_argument('-f', '--bifurcation',
                        action='store_true',
                        help='Print pid if there is a bifurcated sheet.')
    parser.add_argument('-b', '--barrel',
                        action='store_true',
                        help='Print pid if there is a beta barrel.')
    parser.add_argument('-m', '--motif', action='store_true',
                        help='Print metastructure motif.')
    parser.add_argument('-n', '--strands', action='store_true',
                        help='Print number of strands.')
    parser.add_argument('-s', '--save',
                        help='Save results as pkl file in the '
                        'specified directory.')
    parser.add_argument('-d', '--debug', action='store_true')
    args = parser.parse_args()
    main(args.name, args.bifurcation, args.barrel, args.motif, args.strands,
         args.save, args.debug)
