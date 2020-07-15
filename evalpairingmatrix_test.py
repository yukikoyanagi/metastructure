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
    assert math.isclose(epm.compare(pmat, tmat)[0], 1/2)
    assert math.isclose(epm.compare(pmat, tmat)[1], 1/2)
    assert math.isclose(epm.compare(pmat, tmat, upto=2)[0], 1/2)
    assert math.isclose(epm.compare(pmat, tmat, only=1)[0], 0)

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
    
def test_makepartialmatrix():
    cmat = np.array([[ 0, 0,.5,.1, 0],
                     [.8, 0, 0, 0, 0],
                     [ 0,.6, 0,.4, 0],
                     [ 0, 0, 0, 0,.7],
                     [ 0, 0, 0, 0, 0]])
    qmat = np.array([[ 0, 0, 0, 0, 0],
                     [ 1, 0, 0, 0, 0],
                     [ 0, 1, 0, 1, 0],
                     [ 0, 0, 0, 0, 1],
                     [ 0, 0, 0, 0, 0]])
    assert np.allclose(epm.makepmat(cmat), qmat)
    cmat = np.array([[ 0, 0,.2, 0],
                     [ 0, 0, 0,.2],
                     [.1, 0, 0, 0],
                     [ 0,.1,.6, 0]])
    qmat = np.array([[ 0, 0, 1, 0],
                     [ 0, 0, 0, 1],
                     [ 0, 0, 0, 0],
                     [ 0, 0, 1, 0]])
    assert np.allclose(epm.makepmat(cmat), qmat)
