#!/usr/bin/env python3

#   Import required modules from the standard Python library
try: #  Not all versions of Python have argparse, though all versions of Python 3 do
    import argparse
except ImportError:
    sys.exit("Please use Python 3 for this script")

def set_args():
    '''Set the arguments for this example'''
    Arguments = argparse.ArgumentParser(add_help=True, description="A simple Python program that demonstrates brute force algorithms in the context of finding all possible gene sequences given a peptide string", epilog='''To change the peptide string, simply type it after the program name; for example:
        BruteForceRetrotranslate.py MLP''')
    Arguments.add_argument('peptide',
        type=str,
        nargs='?',
        default='FTW',
        help="Enter a peptide string to be used with this example. By default, we use 'FTW'")
    args = Arguments.parse_args()
    return(args)
