#!/usr/bin/env python3

import argparse, random, uuid, pickle

class ReducedSeq(object):

    def __init__(self, seq, length=None):
        '''
        seq is a sequence of CH and CS, possibly ending with a single C
        '''
        if length:
            self.length = length
        else:
            self.length = len(seq)
        if seq[0] == 'C':
            self.seq = self.translate(seq)
        else:
            self.seq = seq

    def __str__(self):
        return self.translate(self.seq)

    @property
    def n_B(self):
        '''
        Number of beta strands
        '''
        return self.seq.count('B')

    @classmethod
    def genseq(cls, length):
        '''
        A random sequence of A's and B's of the given length
        '''
        slength = length//2
        allowed = list(range(slength + 1))
        try:
            del allowed[-2]
        except IndexError: #there is only one allowed value
            pass
        a = random.choice(allowed)
        As = random.sample(range(slength), a)
        seq = ['A' if i in As else 'B' for i in range(slength)]

        return cls(''.join(seq), length=length)

    def translate(self, seq):
        '''
        Translate from (A, B)-sequence to reduced sequence, and vice versa
        '''
        if seq[0] in 'AB':
            newseq = ['CH' if c == 'A' else 'CS' for c in seq]
            if self.length % 2:
                newseq.append('C')
        elif seq[0] == 'C':
            if seq[-1] == 'C':
                seq = seq[:-1]
                newseq = ['A' if c == 'H' else 'B' for c in seq[1::2]]
        return ''.join(newseq)

    def randomgraph(self):
        '''
        Generate random fatgraph based on the sequence
        '''
        while True:
            nums = splitnum(self.n_B)
            if not nums or min(nums) > 1:
                break
        strands = [i for i, c in enumerate(self.seq) if c=='B']
        sheets = []
        for n in nums:
            sheet = random.sample(strands, n)
            sheets.append(sheet)
            strands = [i for i in strands if i not in sheet]

        orientation = [self.randorientation(len(sheet))
                       for sheet in sheets]

        vertices = []
        for sheet, ori in zip(sheets, orientation):
            newsheet = [(i*10, i*10+1) if o > 0 else (i*10+1, i*10)
                        for i, o in zip(sheet, ori)]
            vertices.append(newsheet)

        nvs = []
        for v in vertices:
            if v[0] < v[-1]:
                nvs.append(v)
            else:
                v.reverse()
                nvs.append(v)
        vertices = nvs

        v, e, ie = makefatgraph(vertices)
        
        return v, e, ie

    def randorientation(self, l):
        '''
        Return random sequence of 1 (parallel) and -1 (ant-parallel)
        of length l, starting with 1.
        '''
        return [1] + [random.choice([1, -1]) for i in range(l-1)]

def makefatgraph(vertices):
    '''
    Construct vertices, edges, & internal edges from given vertices.
    '''
    vdict = {}
    for vertex in vertices:
        #We want anti-clockwise ordering of half-edges on each vertex
        l = [v[0] for v in vertex] + [v[1] for v in vertex][::-1]
        r = len(vdict)
        for i in range(r, r+len(l)):
            vdict[l[i - r]] = i+1

    #Now we find edges
    halfedges = sorted([i for vertex in vertices
                        for v in vertex
                        for i in v])
    edges = [(halfedges[i], halfedges[i+1])
             for i in range(1, len(halfedges)-2, 2)]

    #Translate vertices and edges according to vdict
    v = []
    p = 0
    for vertex in vertices:
        v.append(tuple(range(p+1, 2*len(vertex)+p+1)))
        p += 2*len(vertex)

    iv = []
    iedges = [e for vertex in vertices for e in vertex]
    for e in iedges:
        iv.append((vdict[e[0]], vdict[e[1]]))

    e = []
    for d in edges:
        e.append((vdict[d[0]], vdict[d[1]]))

    return v, e, iv

def splitnum(num, lst=[]):
    '''
    Split num into random number of int's, such that the sum = num.
    '''
    if num == 0:
        return lst
    elif num < 4:
        return lst + [num]
    allowed = list(range(2, num+1))
    try:
        del allowed[-2]
    except IndexError:
        pass
    n = random.randint(2, num)
    if n == 0:
        raise ValueError
    return splitnum(num - n, lst + [n])


def main(length, save=None):
    rs = ReducedSeq.genseq(length)
    vertices, edges, in_edges = rs.randomgraph()
    name = str(uuid.uuid4().hex)

    if save:
        with open('{}/{}.pkl'.format(save, name), 'wb') as fh:
            pickle.dump((vertices, edges, in_edges, str(rs)), fh)
    else:
        print('{}:vertices={}'.format(name, vertices))
        print('{}:edges={}'.format(name, edges))
        print('{}:in_edges={}'.format(name, in_edges))
        print('{}:sequence={}'.format(name, str(rs)))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('l', type=int,
                        help='Length of reduced sequence.')
    parser.add_argument('-s', '--save',
                        help='Output directory name for .pkl file')
    args = parser.parse_args()
    main(args.l, args.save)
