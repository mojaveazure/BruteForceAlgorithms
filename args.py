#!/usr/bin/env python3

#   Import required modules from the standard Python library
try: #  Not all versions of Python have argparse, though all versions of Python 3 do
    import argparse
except ImportError:
    sys.exit("Please use Python 3 for this script")

def set_args():
    '''Set the arguments for this example'''
    Arguments = argparse.ArgumentParser(add_help=True)
    Arguments.add_argument('peptide',
        type=str,
        nargs='?',
        default='FTW',
        help="Enter a peptide string to be used with this example")
    args = Arguments.parse_args()
    return(args)
