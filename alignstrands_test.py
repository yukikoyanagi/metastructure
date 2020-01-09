#!/usr/bin/env python3

import numpy as np
import alignstrands as als

def test_makemotif():
    p_mat = np.array([[0., 0.25, 0.35],
                      [0., 0., 0.6],
                      [0., 0., 0.]])
    o_mat = np.array([[False, True, True],
                      [False, False, False],
                      [False, False, False]])
    assert als.makemotif(p_mat, o_mat) == {((1,2), False), ((0,2), True)}
    p_mat = np.array([[0, 1.3, .94, .02, .02, .10, .02],
                      [0, 0,   .37, .02, .02, .05, .02],
                      [0, 0,   0,   .04, .03, .74, .03],
                      [0, 0,   0,   0,   1.9, .04, .02],
                      [0, 0,   0,   0,   0,   .04, .02],
                      [0, 0,   0,   0,   0,   0,   .20],
                      [0, 0,   0,   0,   0,   0,   0]])
    o_mat = np.ones((7, 7), dtype=bool)
    o_mat[3,4] = False
    assert als.makemotif(p_mat, o_mat) == \
        {((0,1), True), ((0,2), True), ((2,5), True),
         ((5,6), True), ((3,4), False), ((1,3), True)}

def test_conn_comp():
    s = {(1,2),(1,3),(3,6),(6,7),(4,5)}
    start = 3
    assert als.conn_comp(s, start) == {(1,2), (1,3), (3,6), (6,7)}
    start = 5
    assert als.conn_comp(s, start) == {(4,5)}
    s.add((2,7))
    print(s)
    start = 3
    assert als.conn_comp(s, start) == {(1,2), (1,3), (3,6), (6,7), (2,7)}

def test_follow():
    s = {(1,2),(1,3),(3,6),(6,7),(4,5)}
    assert als.follow(s, 2) == [2,1,3,6,7]

def test_tovertices():
    mot = {((0,1),True),((0,2),False)}
    assert als.tovertices(mot) == [[(10,11),(0,1),(21,20)],]
    mot = {((1,2), False), ((0,2), False), ((0,5), False), ((3,4), False)}
    assert als.tovertices(mot) == [[(10,11),(21,20),(0,1),(51,50)],[(30,31),(41,40)]]
    mot = {((0,3),True),((0,4),True),((1,2),False)}
    assert als.tovertices(mot) == [[(30,31),(0,1),(40,41)],[(10,11),(21,20)]]
    mot = {((0,3),False),((0,4),False),((1,2),False)}
    assert als.tovertices(mot) == [[(31,30),(0,1),(41,40)],[(10,11),(21,20)]]

def test_hasunpaird():
    mat = np.array([[0,0,0], [0,0,1], [0,0,0]])
    assert als.hasunpaird(mat)
    mat[0,1] = 1
    assert not als.hasunpaird(mat)

def DONOT_test_makematrices():
    strands = [(0,1), (2,3), (4,5)]
    scores = [((0,1),(2,3),0.95,False), ((0,1),(2,3),0.4,False),
              ((0,1),(4,5),0.5,False),  ((0,1),(4,5),0.4,True),
              ((2,3),(4,5),0.6,False)]
    c_mat = np.array([[0,0.95,0.5],[0,0,0.6],[0,0,0]])
    o_mat = np.array([[False, False,  True],
                      [False, False, False],
                      [False, False, False]])
    assert np.array_equal(als.makematrices(strands, scores, 0, 0)[0], c_mat)
    assert np.array_equal(als.makematrices(strands, scores, 1, 0)[1], o_mat)
    c_mat[0,2] = 0
    assert np.array_equal(als.makematrices(strands, scores, 0, 0.6)[0], c_mat)
    c_mat = np.array([[0, 0.4, 0.4],[0,0,0],[0,0,0]])
    assert np.array_equal(als.makematrices(strands, scores, 1, 0)[0], c_mat)
    assert np.array_equal(als.makematrices(strands, scores, 1, 0.5)[0], c_mat)
