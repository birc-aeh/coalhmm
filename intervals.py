def enumerate_all_path_distributions(componentList, noTimePoints):
    '''Creates all sequences of noTimePoints time points assigned
    components componentList in a monotone sequence.'''

    def _sub_seq_generator(pointsLeft, compsLeft):
        if pointsLeft == 0:
            yield ()
        else:
            for i,comp in enumerate(compsLeft):
                for rest in _sub_seq_generator(pointsLeft-1, compsLeft[i:]):
                    yield (comp,) + rest

    return list(_sub_seq_generator(noTimePoints, componentList))

def enumerate_all_path_distributions_intervals(componentList, timePoints):
    '''Similar to the other enumerate, but timePoints is now a list and
    componentList is assumed to be tuples with epochs attached.'''
    nepochs = len(timePoints)
    path_split = [[y for y in componentList if y[0] == i] for i in xrange(nepochs)]
    tmp_paths = [tuple()]
    for i in xrange(nepochs):
        new_paths = []
        for p in tmp_paths:
            extra = enumerate_all_path_distributions(path_split[i], timePoints[i])
            for p2 in extra:
                new_paths.append(p + p2)
        tmp_paths = new_paths
    for p in tmp_paths:
        yield p

def enumerate_all_transitions(SCCPaths, noTimePoints):
    '''Generate all HMM transitions from the list of SCC paths and the
    number of time points.'''
    states = set()
    for path in SCCPaths:
        if isinstance(noTimePoints, list):
            for assignment in enumerate_all_path_distributions_intervals(path,noTimePoints):
                states.add(assignment)
        else:
            for assignment in enumerate_all_path_distributions(path,noTimePoints):
                states.add(assignment)
    return states

if __name__ == "__main__":
    path = [(0,1),(0,2),(1,0),(1,1)]
    nintervals = [3, 2]
    path_split = [[y for y in path if y[0] == i] for i in xrange(len(nintervals))]
    print path_split
    tmp_paths = [tuple()]
    for i in xrange(len(nintervals)):
        new_paths = []
        for p in tmp_paths:
            extra = enumerate_all_path_distributions(path_split[i], nintervals[i])
            for p2 in extra:
                new_paths.append(p + p2)
        tmp_paths = new_paths
    for p in tmp_paths:
        print p

