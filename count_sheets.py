#!/usr/bin/env python3
'''
Count the number of beta sheets for the given protein chain, in the
format '16KPA' (i.e. protein id + chain id). 
Ignore inter-chain bond but still record as a sheet if intra-chain
bonds results in a sheet.
'''

PDB_DIR = '/home/au447708/QGM/metastr/data/pdb'
CHAIN_COLS = [21, 32, 49, 64]
ID_COLS = [11, 12, 13]

import sys, argparse
from collections import Counter

from pdbtools import isinchain

def main(name, height):
    pid = name[:4]
    chain = name[4]

    slst = []

    try:
        with open('{}/{}.pdb'.format(PDB_DIR, pid.lower())) as fh:
            for line in fh:
                if (line.startswith('SHEET') and
                    isinchain(line.strip(), chain)):
                    slst.append(line[ID_COLS[0]:ID_COLS[0] + len(ID_COLS)])
    except FileNotFoundError:
        exit()

    c = Counter(slst)
    sheets = {k:v for k, v in c.items() if v > 1}

    if height:
        if len(sheets) == 0:
            exit()
        print(','.join(map(str,sheets.values())))
    else:
        print(len(sheets))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='Count the number/height of beta-sheets.'
    )
    parser.add_argument('name',
                        help="Protein's name, including chain id.")
    parser.add_argument('-e', '--height',
                        action='store_true',
                        help='Count height of sheets rather than '
                        'number of sheets.')
    args = parser.parse_args()
    main(args.name, args.height)
