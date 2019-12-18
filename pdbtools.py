BIF_FILE = '/home/au447708/QGM/metastr/data/bifurcation.lst'
BAR_FILE = '/home/au447708/QGM/metastr/data/barrels.lst'
PDB_DIR = '/home/au447708/QGM/metastr/data/pdb'

#column numbers in pdb file
CHAIN_COLS = [21, 32, 49, 64]
STRAND_S = 7
STRAND_E = 9
SHEETID_S = 11
SHEETID_E = 13
NUMSTRANDS_S = 14
NUMSTRANDS_E = 15
INITCHAINID = 21
INITSEQNUM_S = 22
INITSEQNUM_E = 25
ENDCHAINID = 32
ENDSEQNUM_S = 33
ENDSEQNUM_E = 36
SENSE_S = 38
SENSE_E = 39

import re
from collections import namedtuple

Strand = namedtuple('Strand', ['num', 'start_c', 'start',
                               'end_c', 'end', 'sense'])

def fromline(line):
    '''
    Construct a Strand namedtuple from a line in pdb file
    '''
    sheet_id = line[SHEETID_S:SHEETID_E+1].strip()
    num = int(line[STRAND_S:STRAND_E+1])
    start_c = line[INITCHAINID]
    end_c = line[ENDCHAINID]
    start = int(line[INITSEQNUM_S:INITSEQNUM_E+1])
    end = int(line[ENDSEQNUM_S:ENDSEQNUM_E+1])
    try:
        sense = int(line[SENSE_S:SENSE_E+1])
    except ValueError: #Some files omit sense when it is 1!!
        sense = 1
    return sheet_id, Strand(num, start_c, start, end_c, end, sense)

def sheetsonchain(name):
    '''
    Return a dict of {sheet_id: [list of Strands]} for the given name.
    Ignore a strand if it is not on the given chain. 
    '''
    pid = name[:4].upper()
    chain = name[4].upper()
    sheets = {}

    with open('{}/{}.pdb'.format(PDB_DIR, pid.lower())) as fh:
        for line in fh:
            if not line.startswith('SHEET'):
                continue
            line = line.strip()
            if not isinchain(line, chain):
                continue
            sheet, strand = fromline(line)
            try:
                sheets[sheet].append(strand)
            except KeyError:
                sheets[sheet] = [strand,]
    return sheets

def isinchain(line, chain):
    '''
    True if the given line of pdb file is in the specified chain
    '''
    cid = [line[i] for i in CHAIN_COLS[:2]]
    return all([(c == chain) for c in cid])

def isregular(pid):
    '''
    False if pid contains bifurcation or beta-barrels.
    '''
    for f in [BIF_FILE, BAR_FILE]:
        with open(f) as fh:
            bad = fh.read()
        if re.search(pid, bad):
            return False
    return True
