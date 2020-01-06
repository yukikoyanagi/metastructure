#!/usr/bin/env python3

import argparse, pickle, pathlib
from itertools import combinations
from math import sqrt

def compare(tmot, pmot, opt):
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
    strands = range(max(s for p in tpairs for s in p) + 1)
    apairs = set(combinations(strands, 2))
    if opt > 0:
        tpairs = {p for p in tpairs if p[1]-p[0]==opt}
        ppairs = {p for p in ppairs if p[1]-p[0]==opt}
        apairs = {p for p in apairs if p[1]-p[0]==opt}
    elif opt < 0:
        opt *= -1
        tpairs = {p for p in tpairs if p[1]-p[0]>opt}
        ppairs = {p for p in ppairs if p[1]-p[0]>opt}
        apairs = {p for p in apairs if p[1]-p[0]>opt}
    TP = len(tpairs & ppairs)
    FP = len(ppairs - tpairs)
    TN = len((apairs - tpairs) & (apairs - ppairs))
    FN = len(tpairs - ppairs)

    c_pairs = set(p for p in pmot if p[0] in tpairs)
    c_links = set(p for p in pmot if p in tmot and p[0] in tpairs)

    return TP, FP, TN, FN, len(c_pairs), len(c_links)

def main(t_dir, p_dir, diag):
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
                TP, FP, TN, FN, CP, CL = compare(tf, pf, diag)
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


    specificity = theTN / (theTN + theFP)
    sensitivity = theTP / (theTP + theFN)
    try:
        correlation = (theTP*theTN - theFP*theFN) / \
                      sqrt((theTP+theFN)*(theTP+theFP)*(theTN+theFN)*(theTN+theFP))
    except ZeroDivisionError:
        correlation = -99
    accuracy = (theTP + theTN) / (theTP + theTN + theFP + theFN)

    print('Specificity: {:.2%}'.format(specificity))
    print('Sensitivity: {:.2%}'.format(sensitivity))
    print('Correlation: {:.4f}'.format(correlation))
    print('Accuracy:    {:.2%}'.format(accuracy))
    try:
        CO = theCL/theCP
    except ZeroDivisionError:
        CO = 0
    print('Correct orientation (for correctly identified pairs): '
          '{:.2%}'.format(CO))


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
    parser.add_argument('-d', '--diag', type=int, default=0,
                        help="Only consider d'th diagonal. If negative "
                        "value is given, consider all entries to the "
                        "right of d'th diagonal.")
    args = parser.parse_args()
    main(args.true_motif_dir, args.pred_motif_dir, args.diag)
