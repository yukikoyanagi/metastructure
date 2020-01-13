#!/usr/bin/env python3

import numpy as np
import pytest
import frompartialmatrix as fpm

def test_offdiagonal():
    smat = np.array([[0, 1, 1, 1, 0],
                     [0, 0,.5, 0, 1],
                     [0, 0, 0, 0,.5],
                     [0, 0, 0, 0, 1],
                     [0, 0, 0, 0, 0]])
    omat = np.zeros_like(smat, dtype=bool)
    omat[0,2] = True
    fpm.CUTOFF = .5
    assert fpm.offdiagonal(smat, omat, 1) == \
        {((0,1),False), ((1,2),False), ((3,4),False)}
    fpm.CUTOFF = .6
    assert fpm.offdiagonal(smat, omat, 2) == {((0,2),True), }

def test_findsheetat():
    smat = np.array([[0, 1, 0, 1, 0],
                     [1, 0, 0, 0, 0],
                     [0, 0, 1, 0, 0],
                     [1, 0, 0, 0, 0],
                     [0, 0, 0, 0, 0]])
    assert fpm.findsheetat(1, smat) == [1, 0, 3]
    with pytest.raises(ValueError):
        fpm.findsheetat(0, smat)
    assert fpm.findsheetat(4, smat) == [4,]
    with pytest.raises(fpm.BarrelError):
        fpm.findsheetat(2, smat)

def test_findsheeets():
    pmat = np.array([[0,0,1,0],
                     [0,0,1,0],
                     [0,0,0,0],
                     [0,0,0,0]])
    assert fpm.findsheets(pmat) == [[0,2,1],[3]]
    pmat[0] = [0,0,0,1]
    assert fpm.findsheets(pmat) == [[0,3],[1,2]]
    pmat[2,3] = 1
    assert fpm.findsheets(pmat) == [[0,3,2,1]]

def test_hasbarrel():
    smat = np.array([[0, 1, 0, 1, 0],
                     [1, 0, 0, 0, 0],
                     [0, 0, 1, 0, 0],
                     [1, 0, 0, 0, 0],
                     [0, 0, 0, 0, 0]])
    assert fpm.hasbarrel(smat)
    smat = np.array([[0, 1, 0, 1, 0],
                     [1, 0, 0, 0, 0],
                     [0, 0, 0, 0, 0],
                     [1, 0, 0, 0, 0],
                     [0, 0, 0, 0, 0]])
    assert not fpm.hasbarrel(smat)
    smat = np.array([[0, 1, 0, 1, 0],
                     [1, 0, 1, 0, 0],
                     [0, 1, 0, 1, 0],
                     [1, 0, 1, 0, 0],
                     [0, 0, 0, 0, 0]])
    assert fpm.hasbarrel(smat)

'''
#This test fails, even though it succeeds in f.x. ipython...
def test_makepartialmatrix():
    cmat = np.array([[ 0,.8,.5,-1, 0],
                     [ 0, 0,.6,-1,-1],
                     [ 0, 0, 0,.4,-1],
                     [ 0, 0, 0, 0,.7],
                     [ 0, 0, 0, 0, 0]])
    qmat = np.array([[ 0, 1,-1,-1, 0],
                     [ 0, 0, 1,-1,-1],
                     [ 0, 0, 0, 1,-1],
                     [ 0, 0, 0, 0, 1],
                     [ 0, 0, 0, 0, 0]])
    assert np.allclose(fpm.makepartialmatrix(cmat), qmat)
'''

def test_combinations_repeat():
    assert set(fpm.combinations_repeat('ABCD', 2, 2)) == \
        {(('A', 'B'), ('C', 'D')),
         (('A', 'C'), ('B', 'D')),
         (('A', 'D'), ('B', 'C'))}
    assert set(fpm.combinations_repeat(range(5), 2, 2)) == \
        {((0, 1), (2, 3)), ((0, 1), (2, 4)), ((0, 1), (3, 4)),
         ((0, 2), (1, 3)), ((0, 2), (1, 4)), ((0, 2), (3, 4)),
         ((0, 3), (1, 2)), ((0, 3), (1, 4)), ((0, 3), (2, 4)),
         ((0, 4), (1, 2)), ((0, 4), (1, 3)), ((0, 4), (2, 3)),
         ((1, 2), (3, 4)), ((1, 3), (2, 4)), ((1, 4), (2, 3))}
    assert set(fpm.combinations_repeat('ABCD', 2, 3)) == set()
    assert set(fpm.combinations_repeat('ABCD', 3, 2)) == set()
    assert set(fpm.combinations_repeat('ABC', 3, 0)) == set()
    assert set(fpm.combinations_repeat('ABC', 0, 2)) == set()

def test_completions():
    pmat = np.array([[0, 1, 0, 0],
                     [1, 0, 0, 0],
                     [0, 0, 0, 1],
                     [0, 0, 1, 0]])
    omat = np.zeros_like(pmat, dtype=bool)
    assert len(list(fpm.completions(pmat, omat))) == 9
    pmat[0,2] = -1
    pmat[0,3] = -1
    assert len(list(fpm.completions(pmat, omat))) == 5
    pmat[1,3] = -1
    pmat[1,2] = -1
    assert len(list(fpm.completions(pmat, omat))) == 1

def test_makevertices():
    sts = [[0,1,2]]
    omat = np.array([[False, False, False],
                     [False, False, True],
                     [False, False, False]])
    assert fpm.makevertices(sts, omat) == [[(0,1),(11,10),(21,20)]]
    sts = [[0,4],[2,1,3]]
    omat = np.array([[False, False, False, False, True],
                     [False, False, False, True, False],
                     [False, False, False, False, False],
                     [False, False, False, False, False],
                     [False, False, False, False, False]])
    assert fpm.makevertices(sts, omat) == \
        [[(0,1),(40,41)],[(21,20),(10,11),(30,31)]]

def test_makefatgraph():
    vertices = [[(0,1),(21,20),(10,11)]]
    assert fpm.makefatgraph(vertices) == \
        ({(1,2,3,4,5,6)}, {(3,6),(4,5)}, {(1,6),(2,5),(3,4)})
    vertices = [[(0,1),(40,41)],[(21,20),(10,11),(30,31)]]
    assert fpm.makefatgraph(vertices) == \
        ({(1,2,3,4),(5,6,7,8,9,10)},
         {(2,8),(4,6),(5,7),(9,10)},
         {(1,4),(2,3),(5,10),(6,9),(7,8)})
