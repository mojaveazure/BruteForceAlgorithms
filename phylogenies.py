#!/usr/bin/env python3

import itertools
from pprint import pprint

#   Looking at solving phylogenetic trees using brute force algorithms
#       First, our sequences
A = 'GATCCATGA'
B = 'GATCTATGC'
C = 'GTCCCATTT'
D = 'AATCCGATC'
E = 'TCTCGATAG'

sample_list = list(A, B, C, D, E)

sampleDict = {
    'A':'GATCCATGA',
    'B':'GATCTATGC',
    'C':'GTCCCATTT',
    'D':'AATCCGATC',
    'E':'TCTCGATAG'
}

def distance_matrix(sampleDict):
    distMat = dict()
    for name2, seq2 in sampleDict.items():
        for name1, seq1 in sampleDict.items():
            if not len(seq1) == len(seq2):
                sys.exit("Whoops!")
            if not name1 == name2:
                count = 0
                for index in range(len(seq2)):
                    if not seq1[index] == seq2[index]:
                        count += 1
                distMat[name1, name2] = count
    return(distMat)


mat = distance_matrix(sampleDict)


def print_distance_matrix(matrix):
    names = list()
    for name1, name2 in matrix.keys():
        names.append(name1)
        names.append(name2)
    names = set(names)
    names = list(names)
    names.sort()
    print(names)
    printMatrix = dict()
    for name1 in names:
        print(name1)
        printMatrix[name1] = []
        print(printMatrix)
        for name2 in names:
            printMatrix[name1].append(matrix.get((name1, name2)))
    for key in printMatrix.keys():
        for index, distance in enumerate(printMatrix[key]):
            printMatrix[key][index] = str(distance)
    print('\t', '\t'.join(names))
    for name in names:
        print(name, '\t', '\t'.join(printMatrix[name]))
    # return(printMatrix)


print_distance_matrix(mat)

[[[['A', 'B'], 'C'], 'D'], 'E']