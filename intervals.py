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

def enumerate_all_transitions(SCCPaths, noTimePoints):
    '''Generate all HMM transitions from the list of SCC paths and the
    number of time points.'''
    states = set()
    for path in SCCPaths:
        for assignment in enumerate_all_path_distributions(path,noTimePoints):
            states.add(assignment)
    return states

if __name__ == "__main__":
    for p in enumerate_all_path_distributions([(0, 1),(0,2),(1,0)], 4):
        print p

