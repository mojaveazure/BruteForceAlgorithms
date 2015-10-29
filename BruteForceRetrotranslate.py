#!/usr/bin/env python3

###################################################
#   Does[0]Compute? 10/29/2015
#   Brute Force Algorithms
#   Paul Hoffman and Diana Trujillo
#   University of Minnesota, Twin Cities
###################################################


#   Import required modules from the standard Python library
import time
import itertools
import re
import sys
import os

#   Import functions from other modules defined in this package
from args import set_args


#   Make sure we're on Python 3
if not sys.version_info[0] == 3:
    sys.exit("Please use Python 3 for this script")

#   Looking at basic brute forcing
#       Predicting gene sequence from protein sequence
#   Dictionary for translating between RNA codons and amino acids
translation = {'GUU': 'V', 'CCG': 'P', 'GCU': 'A', 'CUA': 'L', 'UUA': 'L', 'AUA': 'I', 'ACG': 'T', 'UCG': 'S', 'UUU': 'F', 'GCA': 'A', 'CAG': 'Q', 'GUG': 'V', 'ACU': 'T', 'GGU': 'G', 'GCG': 'A', 'CUU': 'L', 'AAG': 'K', 'UCC': 'S', 'GGA': 'G', 'GGC': 'G', 'CGC': 'R', 'AUC': 'I', 'UUC': 'F', 'ACA': 'T', 'AUG': 'M', 'GGG': 'G', 'CUC': 'L', 'UGU': 'C', 'UCU': 'S', 'AAU': 'N', 'UUG': 'L', 'GAC': 'D', 'GUC': 'V', 'CAC': 'H', 'CAU': 'H', 'GAG': 'E', 'UAU': 'Y', 'CAA': 'Q', 'AGU': 'S', 'UGG': 'W', 'AUU': 'I', 'UAC': 'Y', 'CCA': 'P', 'AGA': 'R', 'CCU': 'P', 'CGU': 'R', 'CGA': 'R', 'AAC': 'N', 'CCC': 'P', 'CGG': 'R', 'CUG': 'L', 'UCA': 'S', 'GCC': 'A', 'GAA': 'E', 'ACC': 'T', 'GAU': 'D', 'AGC': 'S', 'AGG': 'R', 'UGC': 'C', 'AAA': 'K', 'GUA': 'V'}


def check_peptide_sequence(peptide, translation):
    '''Make sure we are given a valid peptide sequence'''
    for amino_acid in peptide: # For each amino acid in our peptide sequence
        if not amino_acid in translation.values(): # If an amino acid in our peptide sequence is invalid
            sys.exit("Invalid amino acid: " + amino_acid) # Exit out with error


def retrotranslate(translation, peptide_length):
    '''Find all possible RNA sequences, given a peptide length'''
    return(''.join(candidate) for candidate in itertools.product(translation.keys(), repeat=peptide_length)) # Create a generator that creates candidate sequences given a peptide length


def RNA_to_peptide(translation, RNA):
    '''Convert RNA sequence into a peptide sequence'''
    peptide = '' # An empty set to hold our complete peptide sequence
    byCodon = re.compile(r'...', re.M) # A regex object to split our RNA sequence into codons
    RNA_codons = byCodon.findall(RNA) # Split our RNA sequence into codons
    for codon in RNA_codons: # For each codon we have
        peptide += translation[codon] # Translate to peptide and append to our existing peptide sequence
    return(peptide) # Return our peptide sequence


def RNA_to_DNA(RNA):
    RNAToCodingDNA = str.maketrans('AUGC', 'ATGC') # Translate table for fining coding DNA sequence from an RNA sequence
    DNA = RNA.translate(RNAToCodingDNA) # Reverse transcribe our RNA sequence to coding DNA
    return(DNA) # Return our reverse transcription


def brute_retrotranslate(translation, peptide):
    '''Brute force DNA sequences that could code for a given peptide'''
    peptide = peptide.upper() # Make sure our peptide sequence is all upper case
    check_peptide_sequence(peptide, translation) # Check to make sure we have a valid peptide sequence
    peptide_length = len(peptide) # How long is our peptide sequence?
    print("Looking at a peptide with " + str(peptide_length) + " amino acids")
    solutions = list() # Create a list to hold our solutions
    print("Generating", len(translation.keys()) ** peptide_length, "candidate sequences...")
    for candidate in retrotranslate(translation, peptide_length): # For each candidate RNA sequence
        generated_peptides = RNA_to_peptide(translation, candidate) # Convert the sequence to peptide
        if generated_peptides == peptide: # If the generated peptide sequence matches
            gene_candidate = RNA_to_DNA(candidate) # Reverse transcribe to the coding DNA sequence
            solutions.append(gene_candidate) # It goes in our solutions list
    print("Found", len(solutions), "possible gene sequences for your peptide")
    return(solutions) # Return our solutions


#       Branch-and-bound version of predicting gene sequence from peptide sequence
def b_and_b_retrotranslate(translation, peptide):
    '''Find candidate DNA sequences using a branch and bound algorithm'''
    peptide = peptide.upper() # Make sure our peptide sequence is all upper case
    check_peptide_sequence(peptide, translation) # Check to make sure we have a valid peptide sequence
    peptide_length = len(peptide) # How long is our peptide sequence?
    print("Looking at a peptide with " + str(peptide_length) + " amino acids")
    solutions = list() # Create a list to hold our solutions
    for i in range(peptide_length): # For each amino acid in the peptide sequence
        holding = list() # Create a list to hold our candidate sequences for this amino acid
        print("Looking at", len(translation.keys()), "possible codons for the current amino acid")
        #for codon in [''.join(candidate) for candidate in itertools.product(translation.keys(), repeat=1)]: # For each possible codon sequence
        for codon in retrotranslate(translation, 1):
            amino_acid = RNA_to_peptide(translation, codon) # Translate the codon to an amino acid
            if amino_acid == peptide[i]: # If we have a match
                gene_codon = RNA_to_DNA(codon) # Reverse translate to the coding DNA sequence
                holding.append(gene_codon) # Add this to our holding list
        print("Generated", len(holding), "candidate codons given this current amino acid")
        if solutions: # If we already have something in our solutions list
            solutions = [''.join(generated) for generated in itertools.product(solutions, holding)] # Add all sequences in our holding list to all sequences in our solutions list
        else: # If our solutions list is empty
            solutions = holding[:] # Copy our holding list to our solutions list
    print("Found", len(solutions), "pssible gene sequences for your peptide")
    return(solutions) # Return our solutions


#   Scoring techniques
def score_retrotranslate(translation, peptide, max_score):
    '''Find candidate DNA sequences by scoring'''
    peptide = peptide.upper() # Make sure our peptide sequence is all upper case
    check_peptide_sequence(peptide, translation) # Check to make sure we have a valid peptide sequence
    peptide_length = len(peptide) # How long is our peptide sequence?
    print("Looking at a peptide with " + str(peptide_length) + " amino acids")
    solutions = list() # Create a list to hold our solutions
    failed = 0 # Let's see how many we failed
    print("Generating", len(translation.keys()) ** peptide_length, "candidate sequences...")
    for candidate in retrotranslate(translation, peptide_length): # For each candidate RNA sequence
        generated_peptides = RNA_to_peptide(translation, candidate) # Convert the sequence to peptide
        score = 0 # Start our score at zero
        for index, amino_acid in enumerate(generated_peptides): # For each amino acid
            if not amino_acid == peptide[index]: # If our generated amino acid does not match our peptide
                score += 1 # Increase our score by one
            if score > max_score: # See if our score exceeds our bounds
                failed += 1 # If so, add one to failed
                break # Break out of our loop
        if score <= max_score:
            gene_candidate = RNA_to_DNA(candidate)
            leadScore = score
            solutions.append(gene_candidate)
    print("We failed", failed, "candidate sequences for your peptide given your scoring bounds")
    print("Found", len(solutions), "possible gene sequences for your peptide")
    return(solutions) # Return our solutions


def brute_time(peptide, translation):
    '''Run our brute force example and time it'''
    print("Searching with a standard brute force algorithm")
    t0 = time.clock() # Start time
    print(brute_retrotranslate(translation, peptide))
    t1 = time.clock() # End time
    bruteTime = t1-t0 # Time difference
    print("The brute force algorithm took", bruteTime, "seconds of walltime")
    return(bruteTime)


def branch_time(peptide, translation):
    '''Run our branch-and-bound example and time it'''
    print("Now running with a branch-and-bound algorithm")
    t0 = time.clock() # Start time
    print(b_and_b_retrotranslate(translation, peptide))
    t1 = time.clock() # End time
    branchTime = t1-t0 # Time difference
    print("The branch-and-bound algorithm took", branchTime, "seconds of walltime")
    return(branchTime)


def score_time(peptide, translation):
    '''Run our scoring example and time it'''
    print("Now running with a scoring algorithm")
    print("We do scoring like golf, high scores are bad")
    print("As the number of mismatches increases, the score increases")
    max_score = int(input("Please enter a maximum scoring bound: ")) # Get a scoring bound and convert it to integer type
    print("Great!")
    t0 = time.clock() # Start time
    print(score_retrotranslate(translation, peptide, max_score))
    t1 = time.clock() # End time
    scoreTime = t1-t0 # Time difference
    print("The scoring algorithm took", scoreTime, "seconds of walltime")
    return(scoreTime)

def clear_screen():
    '''Clear the screen'''
    if os.name == 'nt': # If we're running Windows
        os.system('cls')
    else: # Otherwise, assume Mac/Linux/Unix
        os.system('clear')

def main():
    '''Run the program'''
    args = vars(set_args())
    clear_screen()
    print("Welcome to", sys.argv[0])
    print("The peptide we're using today is", args['peptide'])
    input("Press <enter> to find all gene sequence that could code for this peptide using a brute force algorithm...")
    bruteTime = brute_time(args['peptide'], translation)
    input("Press <enter>")
    repeat = True
    while repeat:
        clear_screen()
        input("Now, let's look at a scoring algorithm, press <enter> to start this...")
        scoreTime = score_time(args['peptide'], translation)
        rep = input("Would you like to repeat the scoring algorithm? If so, type 'yes' without the quotes: ")
        if not rep == 'yes':
            repeat = False
    clear_screen()
    input("Finally, let's look at a branch-and-bound algorithm, press <enter> to start...")
    branchTime = branch_time(args['peptide'], translation)
    input("Press <enter>")
    clear_screen()
    print("To recap:")
    print("Our brute force algorithm took", bruteTime, "seconds of walltime")
    print("Our scoring algorithm took", scoreTime, "seconds of walltime")
    print("Our branch-and-bound algorithm took", branchTime, "seconds of walltime")


main()
