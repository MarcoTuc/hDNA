import networkx
import numpy as np
import pandas as pd

import sys
sys.path.insert(0, '..')
sys.path.insert(0, '.')

import lidakinetics.kineticsfunctions as kf
import lidakinetics.lidakinetics as lk
import lidakinetics.conf as c 

def slide_strands_general(seq1,seq2,min_nucleation):
    l1 = len(seq1)
    l2 = len(seq2)
    n = min_nucleation
    nodes = []
    for b in range(n,l1):
        if b < l1:
            s1 = seq1[0:b]
            s2 = seq2[::-1][0:b][::-1]
            nodes.append([s1,s2])
    nodes.append([seq1,seq2])
    for b in range(1,l1-n+1):
        s1 = seq1[b:]
        s2 = seq2[0:l2-b]
        nodes.append([s1,s2])
    return nodes

def getnodes(strand1,strand2,mincore):
    slids = slide_strands_general(strand1,strand2,mincore)
    nodes = []
    for s in slids:
        a = [lk.iswc(s[0][i],s[1][i]) for i in range(len(s[0]))]
        if core_exists(s,mincore):
            nodes.append(s)
    return nodes


def core_exists(node,coresize):
    nodeA = node[0]
    nodeB = node[1]
    for i in range(len(nodeA) - coresize + 1):
        nA = nodeA[i: i + coresize]
        nB = nodeB[i: i + coresize]
        wcs = sum([lk.iswc(nA[i],nB[i]) for i in range(coresize)])
        if wcs == coresize:
            return True
        break

