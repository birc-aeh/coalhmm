def do_on_all_distributions(S, I, visit):
    buffer = [S[0]] * I
    seen = set()
    def f(S, S_i, I, I_i):
        if S_i == len(S):
            return
        for i in xrange(I_i, I):
            r = buffer[i:]
            buffer[i:] = [S[S_i]] * len(r)
            if buffer not in seen:
                visit(buffer)
                seen.add(buffer)
            f(S, S_i+1, I, i+1)
            buffer[i:] = r
    visit(buffer)
    seen.add(buffer)
    for i in xrange(1,len(S)):
        f(S, i, I, 0)

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
    def p(s):
        print s
    do_on_all_distributions("abc", 4, p)
