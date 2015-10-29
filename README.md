# BruteForceAlgorithms
## A simple Python program for the Does[0]Compute? session on brute force algorithms

This is a simple Python script to demonstrate how brute force algorithms work for [Does\[0\]Compute?](http://morrell-lab.cfans.umn.edu/compute/compute.htm)

We use the task of finding all possible gene sequences that could code for a peptide sequence as our example problem. We find the solution in three ways: a basic brute force algorithm, a brute force algorithm with a scoring mechanism, and a branch-and-bound algorithm.

### Usage
To use, you need [Python 3](https://www.python.org/) installed on your system. If you only have Python 3 on your system, run with the following command:

```shell
python BruteForceRetrotranslate.py
```

If you have Python 2 and Python 3 on your system, run with the following command:
```shell
python3 BruteForceRetrotranslate.py
```

This will use the amino acid sequence of `'FTW'` as the example. To use another sequence, simply type it after the name of the program:

```shell
python BruteForceRetrotranslate.py MLP
```

or

```shell
python3 BruteForceRetrotranslate.py MLP
```
