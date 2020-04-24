#!/usr/bin/env python3

import math
import numpy as np
import evalpairingmatrix as epm

def test_compare():
    pmat = np.array([[0,0,1,0],
                     [0,0,0,0],
                     [0,1,0,0],
                     [0,0,0,0]])
    tmat = np.array([[0,0,1,0],
                     [0,0,1,0],
                     [0,0,0,0],
                     [0,0,0,0]])
    assert math.isclose(epm.compare(pmat, tmat), 5/6)
    assert math.isclose(epm.compare(pmat, tmat, orientation=False), 1)
    assert math.isclose(epm.compare(pmat, tmat, upto=2), 4/5)
    assert math.isclose(epm.compare(pmat, tmat, only=1), 2/3)

def test_makeonematrix():
    pmat = np.array([[0,0,1,0],
                     [0,0,1,0],
                     [0,0,0,0],
                     [0,0,0,0]])
    omat = np.array([[0,0,1,0],
                     [0,0,0,0],
                     [0,0,0,0],
                     [0,0,0,0]])
    omat = omat.astype(bool)
    res = np.array([[0,0,1,0],
                    [0,0,0,0],
                    [0,1,0,0],
                    [0,0,0,0]])
    assert np.allclose(epm.makeonematrix(pmat, omat), res)
    pmat = np.array([[0,0,2,0],
                     [0,0,1,0],
                     [0,0,0,0],
                     [0,0,0,0]])
    omat = np.array([[0,0,0,0],
                     [0,0,0,0],
                     [0,0,0,0],
                     [0,0,0,0]])
    omat = omat.astype(bool)
    res = np.array([[0,0,0,0],
                    [0,0,0,0],
                    [2,1,0,0],
                    [0,0,0,0]])
    assert np.allclose(epm.makeonematrix(pmat, omat), res)

def test_keepupto():
    mat = np.arange(16).reshape((4,4))
    res = np.array([[ 0, 1, 2, 0],
                    [ 4, 5, 6, 7],
                    [ 8, 9,10,11],
                    [ 0,13,14,15]])
    assert np.allclose(epm.keepupto(mat, 2), res)

def test_keeponly():
    mat = np.arange(16).reshape((4,4))
    res = np.array([[ 0, 0, 2, 0],
                    [ 0, 0, 0, 7],
                    [ 8, 0, 0, 0],
                    [ 0,13, 0, 0]])
    assert np.allclose(epm.keeponly(mat, 2), res)
    
