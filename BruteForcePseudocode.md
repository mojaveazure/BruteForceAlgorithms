### Brute Force
```
allCandidates <- generate(Parameters)
solutionSpace <- empty()
for(candidate in allCandidates)
    if valid(candiate, Problem)
        then solutionSpace <- append(candidate)
```

### Reordering Search Space
```
allCandiates <- generate(Parameters)
solutionSpace <- empty()
sortedCandidates <- sort(allCandidates, Problem)
for (candidate in sortedCandidates)
    if valid(candidate, Problem)
        then solutionSpace <- append(candiate)
```

### Branch-and-bound Pseudocode
```
ParameterPart <- split(Parameters)
solutionSpace <- empty()
while not (ParameterPart == null)
    candidates <- generate(ParameterPart)
    for (branch in candidates)
        if not valid(branch, Problem)
            remove(branch)
        solutionSpace <- append(branch)
    ParameterPart <- next(Parameters)
output(solutionSpace)
```

### Scoring Pseudocode
```
allCandidates <- generate(Parameters)
scoringBound <- N
solutionSpace <- empty()
leaderboard <- empty()
for (candidate in allCandidates)
    if valid(candiate, Problem)
        if (score(candidate) > scoringBound)
            leaderboard <- expand(candidate, score(candidate))
leaderboard <- sort(leaderboard by score(candidate))
```

### Help Message
```
usage: BruteForceRetrotranslate.py [-h] [peptide]

A simple Python program that demonstrates brute force algorithms in the
context of finding all possible gene sequences given a peptide string

positional arguments:
  peptide     Enter a peptide string to be used with this example. By default,
              we use 'FTW'

optional arguments:
  -h, --help  show this help message and exit

To change the peptide string, simply type it after the program name; for
example: BruteForceRetrotranslate.py MLP
```