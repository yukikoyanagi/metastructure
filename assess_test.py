#!/usr/bin/env python3

import numpy as np
import assess

def test_findsheeets():
    pmat = np.array([[0,0,1,0],
                     [0,0,1,0],
                     [0,0,0,0],
                     [0,0,0,0]])
    assert assess.findsheets(pmat) == [[0,2,1],[3]]
    pmat[0] = [0,0,0,1]
    assert assess.findsheets(pmat) == [[0,3],[1,2]]
    pmat[2,3] = 1
    assert assess.findsheets(pmat) == [[0,3,2,1]]

def test_makevertices():
    sts = [[0,1,2]]
    omat = np.array([[False, False, False],
                     [False, False, True],
                     [False, False, False]])
    assert assess.makevertices(sts, omat) == [[(0,1),(11,10),(21,20)]]
    sts = [[0,4],[2,1,3]]
    omat = np.array([[False, False, False, False, True],
                     [False, False, False, True, False],
                     [False, False, False, False, False],
                     [False, False, False, False, False],
                     [False, False, False, False, False]])
    assert assess.makevertices(sts, omat) == \
        [[(0,1),(40,41)],[(21,20),(10,11),(30,31)]]

def test_makefatgraph():
    vertices = [[(0,1),(21,20),(10,11)]]
    assert assess.makefatgraph(vertices) == \
        ({(1,2,3,4,5,6)}, {(3,6),(4,5)}, {(1,6),(2,5),(3,4)})
    vertices = [[(0,1),(40,41)],[(21,20),(10,11),(30,31)]]
    assert assess.makefatgraph(vertices) == \
        ({(1,2,3,4),(5,6,7,8,9,10)},
         {(2,8),(4,6),(5,7),(9,10)},
         {(1,4),(2,3),(5,10),(6,9),(7,8)})
