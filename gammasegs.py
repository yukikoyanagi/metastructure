#!/usr/bin/env python3
#
# Extract alpha-gamma segment information for protein
#

#column numbers in summary* file
RESNUM_COL = 1
RES_COL = 2
DSSP_COL = 3
DON_COL = 10
ACC_COL = 11

import argparse, re, pathlib, logging
from itertools import combinations
from operator import itemgetter

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

def getbonds(fname):
    '''
    Return a list of bonds [(from, to),], with each bond expressed
    as a 2-tuple of DSSP residue numbers. 
    '''
    bonds = []
    p = re.compile(r'__[US]{2}')
    with open(fname) as fh:
        for line in fh:
            if p.search(line):
                cols = line.split()
                don = int(cols[DON_COL]) - 1
                acc = int(cols[ACC_COL]) - 1
                bonds.append((don, acc))
    return bonds

def isparallel(f, s, bonds, ext=False):
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
    # Try again with f and s strands extended by 3 residues, but
    # only if they have not been already extended. 
    if ext:
        raise ValueError('Cannot determine parallel/antiparallel.')
    else:
        return isparallel((f[0]-3, f[1]+3), (s[0]-3, s[1]+3),
                          bonds, True)

def isbonded(f, s, bonds):
    '''
    True if strands f and s are bonded
    '''
    for bond in bonds:
        if (bond[0] in range(f[0], f[1]+1) and
            bond[1] in range(s[0], s[1]+1)):
            return True
        if (bond[0] in range(s[0], s[1]+1) and
            bond[1] in range(f[0], f[1]+1)):
            return True
    return False

def main(pid, sd, bd, ar, debug):
    if debug:
        logging.basicConfig(level=logging.DEBUG,
                            format='%(message)s')
    sd = pathlib.Path(sd)
    sf = sd / 'summary{}.txt'.format(pid.upper())
    dseq = dsspseq(sf)
    logging.debug(dseq)
    rseq = resseq(sf)
    logging.debug(rseq)
    bd = pathlib.Path(bd)
    bf = bd / '{}00.txt'.format(pid.upper())
    bonds = getbonds(bf)
    logging.debug(bonds)
    
    strands = []
    for m in re.finditer(r'S+', dseq):
        strands.append((m.start(), m.end()-1))
    logging.debug(strands)
    for strand in strands:
        logging.debug(rseq[strand[0]: strand[1]+1])

    for f, s in zip(strands, strands[1:]):
        if isbonded(f, s, bonds):
            b = True
            try:
                parallel = isparallel(f, s, bonds)
            except ValueError as e:
                # Impossible to determine parallel/antiparallel
                #print('{}({}, {}): {}'.format(pid, f, s, e))
                logging.warning('{}: {} ({},{})'.format(pid, e, f, s))
                continue
        else:
            b = False
            parallel = False
        seg = rseq[f[1]+1 : s[0]]
        if not ar:
            dseg = dseq[f[1]+1 : s[0]]
            lst = []
            for k, s in enumerate(seg):
                if dseg[k] == 'H':
                    lst.append('@')
                else:
                    lst.append(s)
            seg = ''.join(lst)
        print(seg, 'Paired:{}'.format(b), 'Parallel:{}'.format(parallel))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('protid',
                        help='Protein id, including chain id.')
    parser.add_argument('seqdir',
                        help='Directory containing summary* files.')
    parser.add_argument('bonddir',
                        help='Directory containing hbond files.')
    parser.add_argument('-a', '--all-residue', action='store_true',
                        help='Record alpha segment as residues '
                        'instead of as alphas.')
    parser.add_argument('-d', '--debug', action='store_true',
                        help='Verbose output.')
    args = parser.parse_args()
    main(args.protid, args.seqdir, args.bonddir,
         args.all_residue, args.debug)
