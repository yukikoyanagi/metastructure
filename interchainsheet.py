#!/usr/bin/env python3
#
# Check if the given chain contains an inter-chain sheet.
#

PDB_DIR = '/home/au447708/QGM/metastr/data/pdb'

import argparse, re
from pdbtools import sheetsonchain, Strand

def main(name):
    pid = name[:4].upper()
    chain = name[4].upper()
    try:
        sheets = sheetsonchain(name)
    except FileNotFoundError:
        exit()
    for sheet in sheets:
        pat = r'^SHEET\s+[0-9]+\s+{p}\s?(?P<nstr>[0-9]{{1,2}})'.format(p=sheet)
        with open('{}/{}.pdb'.format(PDB_DIR, pid.lower())) as fh:
            pdb = fh.read()
        m = re.search(pat, pdb, re.M)
        try:
            if int(m.group('nstr')) != len(sheets[sheet]):
                print('{}:{}:{}'.format(pid, chain, sheet))
        except AttributeError as e:
            print(name, sheet)
            raise e

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('name',
                        help="Protein's name, including chain id.")
    args = parser.parse_args()
    main(args.name)
