#!/usr/env/bin python

import argparse, pickle
from fatgraph import fatgraph

def main(pkl):
    with open(pkl, 'rb') as fh:
        try:
            g = pickle.load(fh)
        except pickle.UnpicklingError as e:
            print(pkl)
            raise
    try:
        fg = fatgraph.Fatgraph(g[0], g[1])
    except ValueError:
        print(pkl)
        raise
    genus = fg.genus
    if genus < 0:
        raise ValueError('Negative genus for {}.'.format(pkl))
    name = pkl.split('/')[-1].split('.')[0]
    print(len(fg.boundaries), fg.genus, name)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('pklf',
                        help='.pkl file')
    args = parser.parse_args()
    main(args.pklf)
