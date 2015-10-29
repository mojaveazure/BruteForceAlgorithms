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
