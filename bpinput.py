#!/usr/bin/env python3

summary_dir = '/home/au447708/QGM/metastr/data/summaries2020'
fasta_dir = '/home/au447708/QGM/metastr/data/fasta'

import argparse
import pathlib
import re
import json
from operator import itemgetter

def loadseq(seqf):
    seq = ''
    with open(seqf) as fh:
        for line in fh:
            if line[0]=='>':
                continue
            else:
                seq += line
    return seq

def loadfastaseq(pid):
    '''
    Load original pdb sequence from fasta file
    '''
    fdir = pathlib.Path(fasta_dir)
    pid4 = pid.upper()[:-1]
    chain = pid.upper()[-1]
    fs = fdir.glob('{}*'.format(pid4))
    seqf = None

    for f in fs:
        with open(f) as fh:
            for line in fh:
                if chain in line.split('|')[1].split()[1].split(','):
                    seqf = f
                break
        if seqf:
            break
    return loadseq(seqf)

def fill(seq_data, fasta_seq):
    '''
    Fill the gaps in seq_data with the corresponding segments in 
    fasta_seq.
    '''
    seq_data = sorted(seq_data, key=itemgetter(0))
    
    new_seq = []
    for item in seq_data:
        n = item[0] - 1
        l = len(new_seq)
        if n < l:
            new_seq[n] = item
        elif n == l:
            new_seq.append(item)
        else:
            d = n - l
            newseq = [(-1, '.', '.')] * d + [item]
            new_seq += newseq
    while True:
        if new_seq[0][0] < 0:
            _ = new_seq.pop(0)
        else:
            break
    rseq = ''.join([s[1] for s in new_seq])
    sseq = ''.join([s[2] for s in new_seq])
#    print(rseq)
#    print(sseq)
#    print(fasta_seq)

    #Trim 'X' and '.' from start of sequence
    for i, item in enumerate(new_seq):
        if item[1] == 'X' or item[1] == '.':
            continue
        else:
            start = i
            break
    new_seq = new_seq[start:]
    rseq = ''.join([s[1] for s in new_seq])
    sseq = ''.join([s[2] for s in new_seq])
#    print(rseq)

    #Trim fasta_seq
    segs = re.split(r'[.]+|[X]+', rseq)
    m = re.search(segs[0], fasta_seq)
    if m:
        _fasta = fasta_seq[m.start():]
    else:
        #Discrepancy in rseq and fasta. Probably rseq is missing
        #some segment without indicating.
        m = re.search(segs[0][:20], fasta_seq)
        _fasta = fasta_seq[m.start()]
#    print(_fasta)
    
    #Replace  . in residue with corresponding segment in fasta
    m_segs = []
    start = 0
    f_start = 0
    ignore = []
    while True:
        m = re.search(r'[^.X]+(\.+)[^X.]+', rseq[start:])
        if not m:
            break
        pattern = '{}([A-Z]+?){}'.format(
            rseq[start:][m.start():m.start(1)],
            rseq[start:][m.end(1):m.end()])
#        print(pattern)
        start += m.end(1)
        ma = re.search(pattern, _fasta[f_start:])
        if not ma: #'.' in rseq but no match in fasta...
            ignore += list(range(m.start(1), m.end(1)))
            continue
        m_segs.append(ma[1])
        f_start += ma.end(1)

#    print(ignore)
#    print(m_segs)
    flag = 0
    _seq = []
    for i, item in enumerate(new_seq):
        if i in ignore:
            continue
        elif item[1] == '.' and not flag:
            m = m_segs.pop(0)
            for s in m:
                _seq.append((-1,s,'.'))
            flag = 1
        elif item[1] == '.' and flag:
            pass
        else:
            _seq.append(item)
            flag = 0
    new_seq = _seq
    rseq = ''.join([s[1] for s in new_seq])
#    print(rseq)

    #Replace X in residue if a valid letter is available in fasta
    segs = re.split(r'[X.]+', rseq)
    for i, item in enumerate(new_seq):
        if item[1] == 'X':
            new_seq[i] = (item[0], _fasta[i], item[2])
    rseq = ''.join([s[1] for s in new_seq])
#    print(rseq)

    return new_seq

def format_seq(seq):
    '''
    Tidy up secondary structure sequence.
    '''
    sseq = seq.replace('S', 'E')
    newseq = ''
    for i, s in enumerate(sseq):
        if s == '.':
            if sseq[i-1] == 'E' and sseq[i+1] == 'E':
                newseq += 'E'
            elif sseq[i-1] == 'H' and sseq[i+1] == 'H':
                newseq += 'H'
            else:
                newseq += 'C'
        else:
            newseq += s
    return newseq

def main(pid, raw=False):
    sumdir = pathlib.Path(summary_dir)
    fn = sumdir / 'summary{}.txt'.format(pid.upper())

    data = []
    with open(fn) as fh:
        for line in fh:
            cols = line.split()
            resnum = int(cols[1])
            res = cols[2]
            ss = cols[3]
            data.append((resnum, res, ss))

    try:
        fseq = loadfastaseq(pid)
    except TypeError as e:
        print('{}: {}'.format(pid, e))
        raise e

    try:
        seq_data = fill(data, fseq)
    except Exception as e:
        print('{}: {}'.format(pid, e))
        raise e

    if raw:
        outdir = pathlib.Path(raw)
        outf = outdir / '{}.json'.format(pid.upper())
        with open(outf, 'w') as fh:
            json.dump(seq_data, fh)
        exit()

    rseq = ''.join([s[1] for s in seq_data])
    sseq = ''.join([s[2] for s in seq_data])
    sseq = format_seq(sseq)

    print('>{}'.format(pid.upper()))
    print(rseq)
    print(sseq)

if __name__=='__main__':
    parser = argparse.ArgumentParser(
        description='Create betapro input')
    parser.add_argument('pid')
    parser.add_argument('--raw', '-r',
                        help='Save data in raw format (list) in this directory')
    args = parser.parse_args()
    main(args.pid, args.raw)
