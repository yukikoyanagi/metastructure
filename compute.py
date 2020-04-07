#!/usr/bin/env python3

PKL_DIR = '/home/au447708/QGM/metastr/data/fatgraphs_20'

import pathlib, argparse, pickle
from fatgraph import fatgraph
#import fatgraph

def main(name, g, b, s, k, q):
    with open('{}/{}.pkl'.format(PKL_DIR, name), 'rb') as fh:
        vertices, edges, in_edges, seq = pickle.load(fh)
    f = fatgraph.Fatgraph(vertices, edges)
    if g:
        if f.genus < 0:
            raise ValueError('{}. Negative genus.'.format(name))
        else:
            print(len(f.vertices), f.genus, name)
    elif b:
        print(len(f.vertices))
    elif s:
        print('\t'.join([str(len(v)//2) for v in f.vertices]))
    elif k:
        print(sum([len(v) for v in f.vertices])//2, f.genus, name)
    elif q:
        print(name, seq)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('name',
                        help='Protein name, id and chain id.')
    parser.add_argument('-g', '--genus', action='store_true',
                        help='Compute genus for a given protein.')
    parser.add_argument('-b', '--beta', action='store_true',
                        help='Print number of beta sheets.')
    parser.add_argument('-s', '--strands', action='store_true',
                        help='Print number of strands in beta sheets.')
    parser.add_argument('-k', '--chords', action='store_true',
                        help='Print genus and number of half-edges/2.')
    parser.add_argument('-q', '--sequence', action='store_true',
                        help='Print sequence.')
    args = parser.parse_args()
    main(args.name, args.genus, args.beta, args.strands, args.chords,
         args.sequence)
