#!/usr/bin/env python3

import argparse, pickle, pathlib
from itertools import combinations
from math import sqrt

def compare(tmot, pmot):
    with open(tmot, 'rb') as fh:
        thetmot = pickle.load(fh)

    # Motif produced by asfatgraph.py is slightly different format
    # than motif from alignstrands.py. 
    tmot = set()
    for link in thetmot:
        pair, ori = link
        if pair[0] < pair[1]:
            tmot.add(((pair[0], pair[1]), not ori[0]==ori[1]))
        else:
            tmot.add(((pair[1], pair[0]), not ori[0]==ori[1]))
        
    with open(pmot, 'rb') as fh:
        pmot = pickle.load(fh)

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

def main(t_dir, p_dir):
    tdir = pathlib.Path(t_dir)
    pdir = pathlib.Path(p_dir)
    tfs = tdir.glob('*.pkl')
    pfs = pdir.glob('*.pkl')
    theTP = 0
    theFP = 0
    theTN = 0
    theFN = 0
    theCP = 0
    theCL = 0
    for pf in pfs:
        for tf in tfs:
            if pf.name.upper() == tf.name.upper():
                TP, FP, TN, FN, CP, CL = compare(tf, pf)
                theTP += TP
                theFP += FP
                theTN += TN
                theFN += FN
                theCP += CP
                theCL += CL
                break
        else:
            raise FileNotFoundError(
                'Missing motif file: {}'.format(pf.name))


    specificity = theTP / (theTP + theFP)
    sensitivity = theTP / (theTP + theFN)
    correlation = (theTP*theTN - theFP*theFN) / \
                  sqrt((theTP+theFN)*(theTP+theFP)*(theTN+theFN)*(theTN+theFP))

    print('Specificity: {:.2%}'.format(specificity))
    print('Sensitivity: {:.2%}'.format(sensitivity))
    print('Correlation: {:.4f}'.format(correlation))
    print('Correct orientation (for correctly identified pairs): '
          '{:.2%}'.format(theCL/theCP))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Generate various '
                                     'metrics to evaluate pred_motif '
                                     'against true_motif.')
    parser.add_argument('true_motif_dir',
                        help='Directory containing pickle files '
                        'of true motifs.')
    parser.add_argument('pred_motif_dir',
                        help='Directory containing pickle files '
                        'of  predicted motif.')
    args = parser.parse_args()
    main(args.true_motif_dir, args.pred_motif_dir)
