#!/usr/bin/env python3

#   Import required modules from the standard Python library
import time

from BruteForceAlgorthims import check_peptide_sequence
from BruteForceAlgorthims import retrotranslate
from BruteForceAlgorthims import RNA_to_peptide
from BruteForceAlgorthims import RNA_to_peptide

def brute_time(peptide, translation):
    from BruteForceAlgorthims import brute_retrotranslate
    print("Searching with a standard brute force algorithm")
    t0 <- time.clock()
    print(brute_retrotranslate(translation, peptide))
    print(time.clock() - t0)


def branch_time(peptide, translation):
    from BruteForceAlgorthims import b_and_b_retrotranslate
    print("Now running with a branch-and-bound algorithm")
    t0 <- time.clock()
    print(b_and_b_retrotranslate(translation, peptide))
    print(time.clock() - t0)


def score_time(peptide, translation):
    from BruteForceAlgorthims import score_retrotranslate
    