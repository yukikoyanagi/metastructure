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


def test_validorders():
    col = [[0,1,2], [3,4]]
    expected = set([((0,1,2),(3,4)), ((0,2,1),(3,4)), ((1,0,2),(3,4))])
    assert set(mp.validorders(col)) == expected
