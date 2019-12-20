#!/usr/bin/env python3

#column numbers in summary* file
RESNUM_COL = 1
RES_COL = 2
DSSP_COL = 3
DON_COL = 10
ACC_COL = 11

GAP_SCORE = -4
N_ALIGNED = 1000 #choose all matches
MAX_ORDER = 10 # Consider up to n'th score for each strand pair
CUTOFF_VALUES = [k * 0.05 for k in range(6)]
BG_MAT_FILE = 'genbound_arr.pkl'
FACTOR = 8

import argparse, pathlib, re, os, pickle, math
import numpy as np
from collections import Counter
from itertools import filterfalse, combinations, product, chain, islice, groupby
from operator import itemgetter, mul, xor, not_
from functools import reduce

import alignment.substmatrices as substmatrices
import alignment.alignmenthelper as alignmenthelper
import fatgraph.fatgraph as fatgraph

class NoValidMotifError(Exception):
    def __init__(self, pid):
        self.args = ('Cannot generate valid structure for {}. '
                     'Consider raising N_ALIGNED value.'.format(pid))

def prod(iterable):
    return reduce(mul, iterable, 1)

def readrawseq(fname, col):
    '''
    Read col'th column of file fname as a list. 
    '''
    seq = []
    with open(fname) as fh:
        for line in fh:
            cols = line.split()
            n = int(cols[RESNUM_COL])
            c = cols[col]
            if n-1 == len(seq):
                seq.append(c)
            elif n-1 > len(seq):
                seq.extend(['!']*(n-1-len(seq)))
                seq.append(c)
            else:
                seq[n-1] = c
    return seq

def resseq(fname):
    '''
    Get sequence of residues from summary* file. Residue chars are
    ordered by DSSP residue numbering and may contain !. 
    '''
    return ''.join(readrawseq(fname, RES_COL))

def dsspseq(fname):
    '''
    Get sequence of 3-class DSSP from summary* file.
    DSSP chars are ordered by DSSP residue numbering. 
    In addition to H, S, C, the resulting sequence may contain 
    ! for chain break, and/or previous chains (chains B, C, etc. 
    may start at residue number other than 1).
    '''
    seq = readrawseq(fname, DSSP_COL)

    # We need the start residue number
    with open(fname) as fh:
        line = fh.readline()
    cols = line.split()
    start = int(cols[RESNUM_COL]) - 1

    # Single ! sandwiched by two S's is also an S
    for i in range(start, len(seq)):
        if seq[i] == '!' and seq[i-1] == 'S' and seq[i+1] == 'S':
            seq[i] = 'S'
        # ! franked by two S's is also an S
        try:
            if seq[i] == '!' and (
                    (seq[i-1] == 'S' and seq[i-2] != 'S') or
                    (seq[i+1] == 'S' and seq[i+2] != 'S')):
                seq[i] = 'S'
        except IndexError:
            continue

    return ''.join(seq)

def loadlearning(lf):
    lst = []
    with open(lf) as fh:
        for line in fh:
            cols = line.split()
            if cols[-1] == 'True':
                o = True
            else:
                o = False
            lst.append((frozenset(cols[:2]), o))
    return Counter(lst)

def unique_everseen(iterable, key=None):
    "List unique elements, preserving order. Remember all elements ever seen."
    # unique_everseen('AAAABBBCCDAABBB') --> A B C D
    # unique_everseen('ABBCcAD', str.lower) --> A B C D
    seen = set()
    seen_add = seen.add
    if key is None:
        for element in filterfalse(seen.__contains__, iterable):
            seen_add(element)
            yield element
    else:
        for element in iterable:
            k = key(element)
            if k not in seen:
                seen_add(k)
                yield element

def conn_comp(s, start):
    "s is a set of 2-tuples. Return the connected component of s "
    "which contains start."
    # conn_comp({(0,1),(2,3),(1,4),(4,6),(2,5)}, 0) --> {(0,1), (1,4), (4,6)}
    remain = set([i for p in s for i in p])
    comp = set()
    if start in remain:
        comp |= {e for e in s if start in e}
        s = s - comp # Do not mutate; i.e. no '-='
        for i in set([j for p in comp for j in p if j != start]):
            comp |= conn_comp(s, i)
    return comp

def isvalid(motif):
    "True if motif is valid. Checks for barrels and bifurcations"
    links = {s[0] for s in motif}
    # Check for bifurcations
    cnt = Counter(chain.from_iterable(links))
    if cnt.most_common(1)[0][-1] > 2:
        return False
    # Check for barrel
    while cnt.most_common()[-1][-1] < 2:
        start = cnt.most_common()[-1][0]
        links -= conn_comp(links, start)
        if not links:
            return True
        cnt = Counter(chain.from_iterable(links))
    return False

def makemotif(p_mat, o_mat):
    "Construct beta-sheet motif from two matrices"
    w_mat = np.copy(p_mat)
    motif = set()
    while w_mat.max() > 0:
        idx = np.unravel_index(np.argmax(w_mat), w_mat.shape)
        if isvalid(motif | {(idx, o_mat[idx])}):
            motif.add((idx, o_mat[idx]))
        w_mat[idx] = 0

    #Check motif is valid
    strands = set([i for p, o in motif for i in p])
    if len(strands) < len(p_mat):
        raise ValueError('Cannot construct valid motif')
    return motif

def follow(links, start):
    "Follow links from start to end. links is a set of 2-tuples of "
    "integers. Return a list of nodes."
    # follow({(0,1),(1,2),(3,4)},0) --> [0,1,2]
    if not start in chain.from_iterable(links):
        return [start]
    for link in links:
        if start in link:
            thislink = link
            break

    newlinks = links - {thislink}
    newstart = thislink[int(thislink[0]==start)]
    return [start] + follow(newlinks, newstart)

def vertexkey(v):
    return min(chain.from_iterable(v))

def orientvertices(vertices):
    res = []
    for vertex in vertices:
        s = min(vertex)
        if s[0] > s[1]:
            vertex = [(n[1], n[0]) for n in vertex]
        res.append(vertex)
    return res

def tovertices(motif):
    "Translate motif (set of links between beta sheets) to vertices"
    "(list of 2-tuples)"
    # tovertices({((0,1),True),((0,2),False)}) --> [[(10,11),(0,1),(21,20)]]
    themotif = motif.copy()
    links = {s[0] for s in themotif}
    nodes = set(chain.from_iterable(links))
    ends = []
    for node in nodes:
        if len([link for link in links if node in link]) == 1:
            ends.append(node)
    ends.sort()
    vertices = []
    while ends:
        start = ends.pop(0)
        path = follow(links, start)
        ends.remove(path[-1])
        vertex = []
        reverse_next = False
        for node in path:
            vertex.append(tuple(sorted([node*10, node*10+1],
                                       reverse=reverse_next)))
            for elm in themotif:
                if node in elm[0]:
                    reverse_next = not xor(reverse_next, elm[1])
                    themotif.remove(elm)
                    break
        vertices.append(vertex)
    #Fix ordering of vertices
    vertices = sorted(vertices, key=vertexkey)
    #Fix symmetry along vertical for each vertex
    vertices = orientvertices(vertices)
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

def hasunpaird(mat):
    "True if there is i such that i'th row and column are zero"
    for r, c in zip(mat, mat.T):
        if not (r.any() or c.any()):
            return True
    return False

def makematrices(strands, scores, n, m):
    "Make strand pairing matrix and orientation matrix"
    "Use nth score for each pair with cutoff value = m"
    "Cutoff value is only used if it does not result in"
    "unpaired strand."
    c_mat = np.zeros((len(strands), len(strands)))
    o_mat = np.zeros((len(strands), len(strands)), dtype=bool)

    for k, g in groupby(
            sorted(scores, key=itemgetter(0,1), reverse=True),
            key=itemgetter(0,1)):
        i = strands.index(k[0])
        j = strands.index(k[1])
        if j < i:
            i, j = j, i
        try:
            s = sorted(g, key=itemgetter(2), reverse=True)[n]
        except IndexError:
            continue
        c_mat[i, j] = s[2]
        if s[3]:
            o_mat[i, j] = 1

    if hasunpaird(c_mat):
        raise ValueError('Cannot make structure without unpaired strand.')
    
    # Apply cutoff, but do not leave unpaired strand
    entries = []
    for i, j in combinations(range(len(strands)), 2):
        entries.append((i, j, c_mat[i, j]))
    entries = sorted(entries, key=itemgetter(2))
    for entry in entries:
        if entry[2] >= m:
            break
        else:
            row, col = entry[:2]
            w_mat = c_mat.copy()
            w_mat[row, col] = 0
            if hasunpaird(w_mat):
                continue
            else:
                c_mat[row, col] = 0

    return c_mat, o_mat

def adjust(the_motif, score):
    with open(BG_MAT_FILE, 'rb') as fh:
        bg_mat = pickle.load(fh)
    vertices = tovertices(the_motif)
    vs, es, ies = makefatgraph(vertices)
    fg = fatgraph.Fatgraph(vs, es)
    genus = fg.genus
    bound = len(fg.boundaries)
    try:
        score += (math.log1p(bg_mat[genus, bound]) / math.log(466)) * FACTOR
    except IndexError:
        pass
    return score
    
def assign(pids, tid, tcount, sd):
    "Assign pid's to current task, such that tasks have roughly "
    "equal load (sum of square of # of strands)"
    sd = pathlib.Path(sd)
    ps = []
    for pid in pids:
        fn = sd / 'summary{}.txt'.format(pid.upper())
        dseq = dsspseq(fn)
        n = 0
        for k, g in groupby(dseq):
            if k == 'S':
                n += 1
        ps.append((pid, n))
    res = [[] for i in range(tcount)] 
    total = [0 for i in range(tcount)]

    ps = sorted(ps, key=itemgetter(1), reverse=True)
    while ps:
        p = ps.pop(0)
        idx = total.index(min(total))
        res[idx].append(p[0])
        total[idx] += p[1]**2
    return sorted(res[tid], reverse=True)

def compare(tmot, pmot):
    thetmot = tmot

    # Motif produced by asfatgraph.py is slightly different format
    # than motif from alignstrands.py. 
    tmot = set()
    for link in thetmot:
        pair, ori = link
        if pair[0] < pair[1]:
            tmot.add(((pair[0], pair[1]), not ori[0]==ori[1]))
        else:
            tmot.add(((pair[1], pair[0]), not ori[0]==ori[1]))
        
    tpairs = {t[0] for t in tmot}
    ppairs = {p[0] for p in pmot}
    apairs = set(combinations(sorted(set(s for p in tpairs for s in p)),2))
    TP = len(tpairs & ppairs)
    FP = len(ppairs - tpairs)
    TN = len((apairs - tpairs) & (apairs - ppairs))
    FN = len(tpairs - ppairs)

    c_pairs = set(p for p in pmot if p[0] in tpairs)
    c_links = set(p for p in pmot if p in tmot)

    return TP, FP, TN, FN, len(c_pairs), len(c_links)

def align(pid, sd, bd, lf, md, of, mo, debug):
    if debug:
        print('Processing {}.....'.format(pid))
    # Load learning data. llst is a list of unique sequences.
    learn = loadlearning(lf)
    if debug:
        print(learn.most_common(10))
    ll = learn.elements()
    llst = []
    for l in ll:
        llst += list(l[0])
    llst = list(unique_everseen(llst))

    # Load sequences for the given protein id, and find beta strands
    sd = pathlib.Path(sd)
    sf = sd / 'summary{}.txt'.format(pid.upper())
    dseq = dsspseq(sf)
    rseq = resseq(sf)
    if debug:
        print(dseq)
        print(rseq)

    strands = []
    for m in re.finditer(r'S+', dseq):
        strands.append((m.start(), m.end()-1))
    if debug:
        print(strands)
        for strand in strands:
            print(rseq[strand[0]: strand[1]+1])

    # For each strand, compute alignment score for all seqences in llst
    all_scores = {}
    ah = alignmenthelper.AlignmentHelper()
    smat = substmatrices.SubstitutionMatrices()
    for strand in strands:
        this = rseq[strand[0]: strand[1]+1]
        mat = smat.blosum62
        maxscore = ah.alignStrict(
            this.replace('!', '*'), this.replace('!', '*'),
            substMatrix=mat, gapScore=GAP_SCORE, btrace=False)
        scores = []
        for l in llst:
            thisscore = ah.alignStrict(
                this.replace('!', '*'), l.replace('!', '*'),
                substMatrix=mat, gapScore=GAP_SCORE, btrace=False)
            if thisscore > 0:
                scores.append((l, thisscore/maxscore))

        scores = sorted(scores, key=itemgetter(1), reverse=True)
        if debug:
            print(this)
            print(scores[:6])

        all_scores[this] = scores

    # Create a matrix of bond-scores
    learnkeys = {l[0]: l[1] for l in learn.keys()}
    c_scores = []
    for i, j in combinations(strands, 2):
        seq_i = rseq[i[0]: i[1]+1]
        seq_j = rseq[j[0]: j[1]+1]
        top_i = all_scores[seq_i][:N_ALIGNED]
        top_j = all_scores[seq_j][:N_ALIGNED]
        for s, t in product(top_i, top_j):
            key = frozenset([s[0],t[0]])
            if key in learnkeys:
                c_scores.append((i, j, s[1]*t[1], learnkeys[key]))

    if debug:
        print('Strand pairings computed.')

    motifs = []
    for n, m in product(range(MAX_ORDER), CUTOFF_VALUES):
        try:
            c_mat, o_mat = makematrices(strands, c_scores, n, m)
        except ValueError:
            continue
        try:
            motif = makemotif(c_mat, o_mat)
        except ValueError:
            continue
        s_score = sum([c_mat[p[0]] for p in motif])
        motifs.append((motif, s_score))

    if not motifs:
        raise NoValidMotifError(pid)

    if debug:
        print('Candidate motifs generated.')
        print(motifs)

    if md:
        md = pathlib.Path(md)
        with open(md / '{}.pkl'.format(pid), 'rb') as fh:
            tmot = pickle.load(fh)
        thescore = 0
        themot = None
        for mot, _ in motifs:
            TP, FP, TN, FN, NP, NL = compare(tmot, mot)
            denom = math.sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))
            if denom == 0:
                denom = 1
            score = (TP*TN - FP*FN) / denom
            if score > thescore:
                themot = mot
        motif = themot
    else:
        a_motifs = [(mot, adjust(mot, score)) for mot, score in motifs]
        if debug:
            print('Motif scores computed.')
            print(a_motifs)
        motif = max(a_motifs, key=itemgetter(1))[0]
            
    if mo and of:
        of = of.rstrip('/')
        with open('{}/{}.pkl'.format(of, pid), 'wb') as fh:
            pickle.dump(motif, fh)
        return
    if mo:
        print(motif)
        return
        
    vertices = tovertices(motif)
    if debug:
        print(vertices)
    vs, es, ies = makefatgraph(vertices)
    if of:
        of = of.rstrip('/')
        with open('{}/{}.pkl'.format(of, pid), 'wb') as fh:
            pickle.dump((vs, es, ies), fh)
    else:
        print('{}:vertices={}'.format(pid, vs))
        print('{}:edges={}'.format(pid, es))
        print('{}:in_edges={}'.format(pid, ies))

def main(pid, sd, bd, lf, md, of, sl, mo, debug):
    if sl:
        tcount = int(os.environ['SLURM_NTASKS'])
        tid = int(os.environ['SLURM_PROCID'])
        with open(pid) as fh:
            pids = [line.rstrip('\n') for line in fh]
        plst = assign(pids, tid, tcount, sd)
        for pid in plst:
            if not (pathlib.Path(of) / '{}.pkl'.format(pid)).exists():
                try:
                    align(pid, sd, bd, lf, md, of, mo, debug)
                except NoValidMotifError as e:
                    print(e)
    else:
        try:
            align(pid, sd, bd, lf, md, of, mo, debug)
        except NoValidMotifError as e:
            print(e)


if __name__=='__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('pid',
                        help='Protein ID, including chain ID. '
                        'If --slurm is specified it is a file '
                        'containing list of protein id\'s.')
    parser.add_argument('seqdir',
                        help='Directory containing summary* files.')
    parser.add_argument('bonddir',
                        help='Directory containing hbond files.')
    parser.add_argument('learningf',
                        help='File with beta-strand connections.')
    parser.add_argument('-c', '--compare',
                        help='Compare with the true motifs stored '
                        'in the specified directory and choose the '
                        'closest candidate motif.')
    parser.add_argument('-o', '--output',
                        help='Output directory.')
    parser.add_argument('-s', '--slurm', action='store_true',
                        help='Run this on slurm.')
    parser.add_argument('-m', '--motif', action='store_true',
                        help='Return motif instead of fatgraph.')
    parser.add_argument('-d', '--debug', action='store_true',
                        help='Display debug messages.')
    args = parser.parse_args()
    main(args.pid, args.seqdir, args.bonddir, args.learningf,
         args.compare, args.output, args.slurm, args.motif,
         args.debug)
