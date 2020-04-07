#!/usr/bin/env python3

import pickle, argparse, pathlib, random

CUTOFFS = [0.03, 0.025, 0.02, 0.015, 0.01, 0.005]

def main(pdir, root):
    pdir = pathlib.Path(pdir)

    v = CUTOFFS[0]
    fname = '{}_{}.pkl'.format(root, v)
    with open(pdir / fname, 'rb') as fh:
        data = pickle.load(fh)
    remaining = set([d[0] for d in data])
    result = []

    for v in CUTOFFS:
        fname = '{}_{}.pkl'.format(root, v)
        with open(pdir / fname, 'rb') as fh:
            data = pickle.load(fh)
        rdic = dict()
        for item in data:
            pid = item[0] 
            if pid in remaining:
                try:
                    rdic[pid].append(item[1:])
                except KeyError:
                    rdic[pid] = [item[1:]]
        for pid in rdic:
            if sum([p[1] for p in rdic[pid]]) > 0:
                for l in rdic[pid]:
                    result.append((pid, *l))
                remaining.remove(pid)
        print('Remaining after {}: {}'.format(v, remaining))

    if remaining:
        # Proteins without accepted structure. Choose one random.
        for pid in remaining:
            l = len(rdic[pid])
            x = random.choice(range(l))
            for i in range(l):
                item = rdic[pid][i]
                if i==x:
                    result.append((pid, item[0], True, *item[2:]))
                else:
                    result.append((pid, *item))

    fname = '{}_all.pkl'.format(root)
    with open(pdir / fname, 'wb') as fh:
        pickle.dump(result, fh)

if __name__=='__main__':
    parser = argparse.ArgumentParser(
        description='Select raw results from different cutoff values '
        'so that all proteins have at least one accepted candidate. '
        'The cutoff values are {}.'.format(CUTOFFS)
    )
    parser.add_argument('pkl_dir',
                        help='Directory for pickle files.')
    parser.add_argument('fname_root',
                        help='Filename root. The script will look at '
                        'all files of the form FNAME_ROOT_cutoff.pkl')
    args = parser.parse_args()
    main(args.pkl_dir, args.fname_root)
