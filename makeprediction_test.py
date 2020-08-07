#!/usr/bin/env python3

import numpy as np
import pytest
import makeprediction as mp

def test_allpmats():
    n = 3
    expected = list()
    expected.append(np.array([[0,0,0],[1,0,0],[0,1,0]]))
    expected.append(np.array([[0,0,0],[1,0,1],[0,0,0]]))
    expected.append(np.array([[0,1,0],[0,0,0],[0,1,0]]))
    expected.append(np.array([[0,1,0],[0,0,1],[0,0,0]]))
    expected.append(np.array([[0,0,0],[0,0,0],[1,1,0]]))
    expected.append(np.array([[0,0,0],[0,0,1],[1,0,0]]))
    expected.append(np.array([[0,0,1],[0,0,0],[0,1,0]]))
    expected.append(np.array([[0,0,1],[0,0,1],[0,0,0]]))
    expected.append(np.array([[0,0,0],[1,0,0],[1,0,0]]))
    expected.append(np.array([[0,0,1],[1,0,0],[0,0,0]]))
    expected.append(np.array([[0,1,0],[0,0,0],[1,0,0]]))
    expected.append(np.array([[0,1,1],[0,0,0],[0,0,0]]))

    pmats = list(mp.allpmats(n))
    assert len(pmats) == 12

    assert len(list(mp.allpmats(4))) == 108
    assert len(list(mp.allpmats(4, paired=[(3,2),]))) == 26
    assert len(list(mp.allpmats(4, paired=[(3,2),(1,0)]))) == 9

def test_validorders():
    col = [[0,1,2], [3,4]]
    expected = set([((0,1,2),(3,4)), ((0,2,1),(3,4)), ((1,0,2),(3,4))])
    assert set(mp.validorders(col)) == expected

def test_select_pairs():
    mat = np.array([[0.0, 1.8, 0.0, 2.2],
                    [0.0, 0.0, 0.3, 0.0],
                    [1.3, 0.8, 0.0, 0.1],
                    [2.0, 1.5, 0.5, 0.0]])
    assert mp.select_pairs(mat, 30) == [(0,3),(0,1),(2,1)]
    mat = np.array([[0.0, 2.8, 0.0, 0.2, 0.0, 0.0],
                    [0.0, 0.0, 1.9, 0.0, 0.0, 0.0],
                    [1.3, 0.8, 0.0, 2.7, 0.0, 0.0],
                    [1.0, 0.5, 1.5, 0.0, 0.0, 0.0],
                    [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                    [0.0, 0.0, 0.0, 0.0, 0.3, 0.0]])
    assert mp.select_pairs(mat, 15) == [(0,1),(2,3),(1,2),(5,4)]

def test_pairs2segs():
    ps = [(0,1),(1,2),(3,4)]
    assert mp.pairs2segs(ps) == [[0,1,2],[3,4]]
    ps = [(1,0),(2,1),(4,3)]
    assert mp.pairs2segs(ps) == [[0,1,2],[3,4]]
    ps = [(0,1),(2,3),(4,5),(3,0)]
    assert mp.pairs2segs(ps) == [[1,0,3,2],[4,5]]
