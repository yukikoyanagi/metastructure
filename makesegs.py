#!/usr/bin/env python3

rawseq_dir = '/home/au447708/QGM/metastr/data/rawseq'
bond_dir = '/home/au447708/QGM/metastr/data/HQ60_20'

DON_COL = 10
ACC_COL = 11

import argparse
import pathlib
import json
import re
import logging
from operator import itemgetter

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

def fixdssp(seq):
    '''
    Fix DSSP sequence
    '''
    _seq = seq.replace('S', 'E')
    newseq = ''
    for i, s in enumerate(_seq):
        if s == '.':
            if _seq[i-1] == 'E' and _seq[i+1] == 'E':
                newseq += 'E'
            elif _seq[i-1] == 'H' and _seq[i+1] == 'H':
                newseq += 'H'
            else:
                newseq += 'C'
        else:
            newseq += s
    return newseq

def replaceseq(data, f, seq):
    if len(data) != len(seq):
        raise ValueError('Lengths of two sequences do not match.')
    _data = []
    for d, s in zip(data, seq):
        _d = list(d)
        _d[f] = s
        _data.append(tuple(_d))
    return _data

def isbonded(f, s,  bonds):
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

def getresnums(data, k):
    return (data[k[0]][0], data[k[1]][0])

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
    # If it still fails, we say it's parallel.
    if ext:
        return True
    else:
        return isparallel((f[0]-3, f[1]+3), (s[0]-3, s[1]+3),
                          bonds, True)


def main(pid, skip, debug):
    if debug:
        logging.basicConfig(level=logging.DEBUG,
                            format='%(message)s')

    bd = pathlib.Path(bond_dir)
    bf = bd / '{}00.txt'.format(pid.upper())
    bonds = getbonds(bf)
    logging.debug(bonds)

    rd = pathlib.Path(rawseq_dir)
    rf = rd / '{}.json'.format(pid.upper())
    with open(rf, 'r') as fh:
        seq_data = json.load(fh)
    logging.debug(seq_data)

    dseq = ''.join([s[2] for s in seq_data])
    dseq = fixdssp(dseq)
    logging.debug(dseq)
    seq_data = replaceseq(seq_data, 2, dseq)
    rseq = ''.join([s[1] for s in seq_data])

    #Get strands (as indices for seq_data)
    strands = []
    for m in re.finditer(r'E+', dseq):
        strands.append((m.start(), m.end()-1))
    logging.debug(strands)

    for first, second in zip(strands, strands[skip+1:]):
        #We need residue numbers
        f = getresnums(seq_data, first)
        s = getresnums(seq_data, second)
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

        seg = rseq[first[1]+1 : second[0]]
        dseg = dseq[first[1]+1 : second[0]]
        lst = []
        for k, s in enumerate(seg):
            if dseg[k] == 'H':
                lst.append('@')
            elif dseg[k] == 'E':
                lst.append('=')
            else:
                lst.append(s)
        seg = ''.join(lst)
        print('0', seg, 'Paired:{}'.format(b), 'Parallel:{}'.format(parallel))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('pid')
    parser.add_argument('-s', '--skip', type=int, default=0,
                        help='Include s beta strands in the resulting '
                        'segments. Beta strands are shown as =.')
    parser.add_argument('-d', '--debug', action='store_true',
                        help='Verbose output')
    args = parser.parse_args()
    main(args.pid, args.skip, args.debug)
