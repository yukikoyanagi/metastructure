#!/usr/bin/env python3
#

import argparse, json, logging, pickle
from itertools import permutations, product
from collections import Counter
from collections.abc import Iterable
import numpy as np

from fatgraph.fatgraph import FatgraphB

MAT_FILE = 'gbs_arr_20.pkl'

class AddLinkError(ValueError):
    pass

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

def pairsapart(seq, pairs):
    "True if any of the pairs are apart in the seq"
    for pair in pairs:
        i, j = pair
        if abs(seq.index(i)-seq.index(j))>1:
            return True
    return False

def permutations2(iterable, r=None):
    "Return normal itertools.permutations, if iterable does not "
    "contain list. If it does, keep the elements in list together."
    '''
    [0,1,2] --> 012 021 102 120 201 210
    [0,[1,2]] --> 012 021 120 210
    '''
    flag = False
    for i in iterable:
        if isinstance(i, Iterable):
            flag = True
            break

    if flag:
        for p in permutations(iterable):
            segs = [isinstance(x, Iterable) for x in p]
            for seq in product([1, -1], repeat=sum(segs)):
                out = []
                for k, t in enumerate(segs):
                    if t:
                        i = seq[sum(segs[:k])]
                        out.extend([x for x in p[k][::i]])
                    else:
                        out.append(p[k])
                yield tuple(out)
    else:
        for p in permutations(iterable):
            yield p

def validorders(collection):
    "Return all valid reorderings of elements in collection. "
    "Collection is a list of lists (may contain segments as lists)."
    " Rules for valid orders are; "
    "-within each list, the first element is smaller than the last"
    "-lists are ordered by the smallest element in each list"
    ls = []
    for c in collection:
        l = []
        for p in permutations2(c):
            if p[0] < p[-1]:
                l.append(p)
        ls.append(l)
    ls = sorted(ls, key=min)
    return product(*ls)

def findpidx(ps, ks):
    "Return list of indices in ps where element from ks is found"
    ids = []
    for k in ks:
        try:
            i = ps.index(k)
        except ValueError:
            i = ps.index(k[::-1])
        ids.append(i)

    return ids

def preset_ornt(ornts, ids, ps):
    "Return true if ornts at ids match those in ps"
    for k in range(len(ids)):
        i = ids[k]
        if ornts[i] != (ps[k][0]<ps[k][1]):
            return False
    return True

def pairs2segs(pairs):
    "Convert list of pairs to list og segments"
    '''
    [(0,1),(1,2),(3,4)] --> [[0,1,2],[3,4]]
    '''
    segs = []
    for p in pairs:
        i, j = p
        matched = [s for s in segs if i in s or j in s]
        if len(matched) == 0:
            segs.append([min(p), max(p)])
        elif len(matched) == 1:
            s = matched[0]
            if i == s[0]:
                news = [j, i] + s[1:]
            elif i == s[-1]:
                news = s[:-1] + [i, j]
            elif j == s[0]:
                news = [i, j] + s[1:]
            elif j == s[-1]:
                news = s[:-1] + [j, i]
            segs.remove(s)
            segs.append(news)
        elif len(matched) == 2:
            s, t = matched
            if set([s[0], t[0]]) == set([i,j]):
                news = s[::-1] + t
            elif set([s[0], t[-1]]) == set([i,j]):
                news = s[::-1] + t[::-1]
            elif set([s[-1], t[0]]) == set([i,j]):
                news = s + t
            elif set([s[-1], t[-1]]) == set([i,j]):
                news = s + t[::-1]
            segs.remove(s)
            segs.remove(t)
            segs.append(news)

    for s in segs:
        if s[0] > s[-1]:
            s.reverse()
    
    return sorted(segs, key=min)

def allpmats(n, paired=[]):
    "Yield all valid, complete pairing matrices of size n"
    "paired is a list of tuples representing paired strands. "
    "For a pair (i,j), if i<j the pairing is parallel, if i>j "
    "the pairing is anti-parallel."
    collection = list(range(n))
    if paired:
        segs = pairs2segs(paired)
        for seg in segs:
            for i in seg:
                collection.remove(i)
            collection.append(seg)

    logging.debug('Sheets: {}'.format(collection))

    for p in partition(collection):
        if paired:
            #If paired was specified, p contains pairs as list
            newp = []
            for c in p:
                newc = []
                for x in c:
                    if isinstance(x, Iterable):
                        newc += [y for y in x]
                    else:
                        newc.append(x)
                newp.append(newc)
        else:
            newp = p
        if min([len(s) for s in newp]) < 2:
            continue
        y = 0
        for q in validorders(p):
            y += 1
            if y % 100 == 0:
                logging.debug('Processing {} in {}'.format(y, p))
            #q is a tuple of tuples representing strand pairings
            links = [(i,j) for r in q for i,j in zip(r, r[1:])]
            p_idx = findpidx(links, paired)
            for oris in product([True, False],
                                     repeat=len(links)-len(paired)):
                oris = list(oris)
                link_oris = []
                for k, link in enumerate(links):
                    if k in p_idx:
                        pair = paired[p_idx.index(k)]
                        link_oris.append(pair[0]<pair[1])
                    else:
                        link_oris.append(oris.pop(0))
                mat = np.zeros((n, n))
                for idx, par in zip(links, link_oris):
                    if par:
                        mat[min(idx), max(idx)] = 1
                    else:
                        mat[max(idx), min(idx)] = 1
                yield mat

def addlink(sheets, link):
    "Return sheets (list of lists, each sheet may contain multiple "
    "occurrances of a strand) with the new link added."
    "Gives AddLinkError if adding the link results in barrel or "
    "bifurcation."
    i, j = link
    added = False
    for sheet in sheets:
        if i in sheet or j in sheet:
            sheet.extend([i,j])
            cnt = Counter(sheet).values()
            if list(cnt).count(1) != 2 or list(cnt).count(3) > 0:
                #either barrel (<2) or bifurcation (>2)
                del sheet[-2:]
                raise AddLinkError
            added = True
            break

    if added:
        #check sheets
        if len(sheets) > 1:
            cnt = Counter([x for s in sheets for x in set(s)])
            k, c = cnt.most_common(1)[0]
            if c > 1: #two of the sheets should be merged
                a, b = [s for s in sheets if k in s]
                newsheet = a + b
                #Check newsheet for barrel or bifurcation
                cnt = Counter(newsheet).values()
                if list(cnt).count(1) != 2 or list(cnt).count(3) > 0:
                    if i in a and j in a:
                        a.remove(i)
                        a.remove(j)
                    else:
                        b.remove(i)
                        b.remove(j)
                    raise AddLinkError
                else:
                    sheets.remove(a)
                    sheets.remove(b)
                    sheets.append(newsheet)
    else:
        #No sheet contains i nor j
        sheets.append([i,j])
    return sheets

def select_pairs(mat, p):
    "Select top p% entries from mat and retun as a list of tuples, "
    "while avoiding bifurcation, barrel."
    wmat = np.copy(mat)
    k = wmat.shape[0]
    n = int((k**2 - k) * p/100)
    if n > k-1:
        n = k-1
    sheets = []
    pairs = []
    while n > len(pairs):
        i,j = np.unravel_index(np.argmax(wmat), wmat.shape)
        pairs.append((i,j))
        try:
            sheets = addlink(sheets, (i,j))
        except AddLinkError:
            pairs.remove((i,j))
#        logging.debug(sheets)
#        logging.debug(pairs)
        wmat[i,j] = 0
    return pairs

def main(pf, save, alpha, top, dbg):
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
    paired = select_pairs(pmat, top)

    logging.debug(paired)

    sums = {}
    for cmat in allpmats(n, paired):
        if beta > 0:
            fg = FatgraphB.from_pmat(cmat)
            g = fg.genus
            b = len(fg.boundaries)
            k = max(len(v) for v in fg.vertices)//2
            try:
                s = sums[k]
            except KeyError:
                s = np.sum(bg_mat[:,:,k])
                sums[k] = s
            if s == 0:
                topo_score = 0
            else:
                try:
                    topo_score = bg_mat[g,b,k] / s
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
    parser.add_argument('-t', '--top', type=float, default=0.0,
                        help='Use top t%% of entries from the input '
                        'pairing matrix.')
    parser.add_argument('-d', '--debug', action='store_true',
                        help='Print debug messages.')
    args = parser.parse_args()
    main(args.pmat_file, args.save, args.alpha, args.top, args.debug)
