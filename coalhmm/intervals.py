from itertools import combinations_with_replacement, groupby
def enumerate_all_path_distributions_intervals(componentList, timePoints):
    '''Similar to the other enumerate, but timePoints is now a list and
    componentList is assumed to be tuples with epochs attached.'''
    nepochs = len(timePoints)
    path_split = [list(g) for k,g in groupby(componentList, key=lambda v: v[0])]
    tmp_paths = [tuple()]
    for i in xrange(nepochs):
        new_paths = []
        for p in tmp_paths:
            next_p = path_split[i]
            first = next_p[0]
            extra = combinations_with_replacement(next_p, timePoints[i]-1)
            for p2 in extra:
                new_paths.append(p + (first,) + p2)
        tmp_paths = new_paths
    for p in tmp_paths:
        yield p

def enumerate_all_transitions(SCCPaths, noTimePoints):
    '''Generate all HMM transitions from the list of SCC paths and the
    number of time points.'''
    assert isinstance(noTimePoints, list)
    states = set()
    for path in SCCPaths:
        for assignment in enumerate_all_path_distributions_intervals(path,noTimePoints):
            states.add(assignment)
    return sorted(states)
