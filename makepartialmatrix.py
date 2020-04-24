#!/usr/bin/env python3

#column numbers in summary* file
RESNUM_COL = 1
RES_COL = 2
DSSP_COL = 3
DON_COL = 10
ACC_COL = 11

GAP_SCORE = -4
N_ALIGNED = 10 #choose all matches
MAX_ORDER = 1 # Consider up to n'th score for each strand pair

MATCH = 4
MISMATCH = -4

import argparse, pathlib, re, os, pickle, logging
import multiprocessing as mp
import numpy as np
from datetime import datetime
from itertools import groupby
from operator import itemgetter

import alignment.substmatrices as substmatrices
import alignment.alignmenthelper as alignmenthelper

# Global variables
learn = []
sumdir = pathlib.Path('.')
cpus = None

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

def getd(seg):
    "Return number of sheets ('=' segment)in seg"
    return len([s for s in seg.split('=') if s]) - 1

def loadsegfile(sf, n):
    def istrue(s):
        if s.split(':')[1] == 'True':
            return True
        else:
            return False

    sd = {}
    with open(sf) as fh:
        for line in fh:
            cols = line.split()
            s = cols[1]
            v = [istrue(cols[i]) for i in [2, 3]]
            try:
                [sd[s][i].append(v[i]) for i in range(2)]
            except KeyError:
                sd[s] = [[v[0],], [v[1],]]
    segdicts = [dict() for i in range(n+1)]
    for k in sd:
        d = getd(k)
        if d > n:
            continue
        b = sum(sd[k][0]) / len(sd[k][0])
        if sum(sd[k][0]) == 0:
            o = False
        else:
            o = (sum(sd[k][1]) / sum(sd[k][0]) >= 0.5)
        segdicts[d][k] = (b, o)
    return segdicts

def getstrands(seq):
    strands = []
    for m in re.finditer(r'S+', seq):
        strands.append((m.start(), m.end()-1))
    return strands

def togammaseq(r_seq, d_seq):
    "Given sequences of residues (r_seq) and dssp classes (d_seq), "
    "return a sequence where helices are '@' and strands are '=', "
    "and the rest is shown as residues."
    seq = []
    for k, s in enumerate(r_seq):
        if d_seq[k] == 'H':
            seq.append('@')
        elif d_seq[k] == 'S':
            seq.append('=')
        else:
            seq.append(s)
    return ''.join(seq)

def makematrices(scores):
    "Make pairing and orientation matrices from scores, which is a "
    "list of tuples (pid, segment, [candidates], alignment_score)."
    "Set pairing score to -1 if the corresponding strands should "
    "not be paired."

    logger = logging.getLogger()
    logger.debug('Making matrices from:\n{}'.format(scores))
    pid = scores[0][0]
    sf = sumdir / 'summary{}.txt'.format(pid.upper()) #sumdir is global
    dseq = dsspseq(sf)
    rseq = resseq(sf)
    aseq = togammaseq(rseq, dseq)
    strands = getstrands(dseq)
    s_mat = np.zeros((len(strands), len(strands)))
    o_mat = np.zeros((len(strands), len(strands)), dtype=bool)
    for score in scores:
        pid, seg, cands, a_score = score
        k = getd(seg)
        try:
            p_score = sum([learn[k][c][0] for c in cands]) / len(cands)
        except ZeroDivisionError:
            p_score = 0
        s_score = p_score * a_score
        if s_score == 0:
            s_score = -1
        try:
            o_score = sum([learn[k][c][1] for c in cands]) / len(cands)
        except ZeroDivisionError:
            o_score = 0
        o = o_score >= 0.5
        for s, t in zip(strands, strands[k+1:]):
            if seg == aseq[s[1]+1: t[0]]:
                i = strands.index(s)
                j = strands.index(t)
                s_mat[i,j] = s_score
                o_mat[i,j] = o
    logger.debug('Constructed matrices for {}.'.format(pid))
    return s_mat, o_mat

def align(seg, candidates, pid):
    'Find seq(s) in candidates that best aligns with seg.'
    'Return the best candidate and normalised alignment score.'
    logger = logging.getLogger()
    start = datetime.now()
    logger.debug('{}: Aligning {}'.format(start, seg))
    ah = alignmenthelper.AlignmentHelper()
    smat = substmatrices.SubstitutionMatrices()
    mat = smat.blosum62plus2(match=MATCH, mismatch=MISMATCH)
    maxscore = ah.alignStrict(
        seg.replace('!', '*'), seg.replace('!', '*'),
        substMatrix=mat, gapScore=GAP_SCORE, btrace=False
    )
    currentmax = 0
    bestcands = []
    for cand in candidates:
        if cand == seg:
            currentmax = maxscore
            bestcands = [cand,]
            break
        m = max(mat.values())
        if min(len(seg),len(cand))*m - \
           abs(len(seg)-len(cand))*GAP_SCORE < 0:
            continue
        score = ah.alignStrict(
            seg.replace('!', '*'), cand.replace('!', '*'),
            substMatrix=mat, gapScore=GAP_SCORE, btrace=False
        )
        if score < currentmax:
            continue
        elif score > currentmax:
            bestcands = [cand,]
            currentmax = score
        else:
            bestcands.append(cand)
    end = datetime.now()
    d = end - start
    logger.debug('{}: Finished {}. Elapsed {}'.format(end, seg, d.seconds))
    return pid, seg, bestcands, currentmax/maxscore

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
    return res[tid]

def make(pids, sd, lf, of, up):
    logger = logging.getLogger()
    logger.debug('Processing {}....'.format(';'.join(pids)))

    global learn
    learn = loadsegfile(lf, up)
    global sumdir
    sumdir = pathlib.Path(sd)

    segs = []
    for pid in pids:
        logger.debug('Loading {}....'.format(pid))
        sf = sumdir / 'summary{}.txt'.format(pid.upper())
        dseq = dsspseq(sf)
        rseq = resseq(sf)
        logger.debug(dseq)
        logger.debug(rseq)

        strands = getstrands(dseq)
        logger.debug(strands)

        for i in range(up):
            for s, t in zip(strands, strands[i+1:]):
                dseg = dseq[s[1]+1 : t[0]]
                rseg = rseq[s[1]+1 : t[0]]
                seg = togammaseq(rseg, dseg)
                segs.append((seg, i, pid))

    logger.debug(segs)
    logger.debug('Aligning {} segments...'.format(len(segs)))

    segs = sorted(segs, key=lambda s: len(s[0]), reverse=True)

    pool = mp.Pool(processes=cpus)
    res_objects = [pool.apply_async(
        align, args=(seg, list(learn[i].keys()), pid))
        for seg, i, pid in segs]
    pool.close()
    pool.join()
    scores = [r.get() for r in res_objects]

    logger.debug(scores)
    logger.debug('Obtained {} scores.'.format(len(scores)))

    for pid, g in groupby(
            sorted(scores, key=itemgetter(0)),
            key=itemgetter(0)):
        s_mat, o_mat = makematrices(list(g))
        logger.debug(pid)
        logger.debug(s_mat)
        logger.debug(o_mat)
        if of:
            ofile = pathlib.Path(of) / '{}.mat.pkl'.format(pid)
            with open(ofile, 'wb') as fh:
                pickle.dump((s_mat, o_mat), fh)
        else:
            print(pid)
            print(s_mat)
            print(o_mat)

def main(pids, sd, lf, of, sl, up, pr, debug):
    if debug:
        logging.basicConfig(format='%(message)s', level=logging.DEBUG)
    if of:
        #Check the output directory exists
        od = pathlib.Path(of)
        if not od.exists():
            raise FileNotFoundError('{} not found.'.format(of))
    if pr:
        global cpus
        cpus = pr
    if pathlib.Path(pids).is_file():
        with open(pids) as fh:
            pids = [line.rstrip('\n') for line in fh]
    else:
        pids = pids.split(',')
    if sl:
        tcount = int(os.environ['SLURM_NTASKS'])
        tid = int(os.environ['SLURM_PROCID'])
        plst = assign(pids, tid, tcount, sd)
        make(plst, sd, lf, of, up)
    else:
        make(pids, sd, lf, of, up)


if __name__=='__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('pids',
                        help='Comma-separated list of protein IDs, '
                        'including chain ID. Or a file '
                        'containing list of protein id\'s.')
    parser.add_argument('seqdir',
                        help='Directory containing summary* files.')
    parser.add_argument('learningf',
                        help='File with alpha/gamma segments.')
    parser.add_argument('-o', '--output',
                        help='Output directory.')
    parser.add_argument('-s', '--slurm', action='store_true',
                        help='Run this on slurm.')
    parser.add_argument('-u', '--upto', type=int,
                        help='Construct pairing matrix up to u\'th '
                        'diagonal.')
    parser.add_argument('-p', '--processes', type=int,
                        help='Only use p processes to compute alignment.')
    parser.add_argument('-d', '--debug', action='store_true',
                        help='Display debug messages.')
    args = parser.parse_args()
    main(args.pids, args.seqdir, args.learningf,
         args.output, args.slurm, args.upto, args.processes, args.debug)
