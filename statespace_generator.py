
from sets import ImmutableSet as iset

class CoalSystem(object):
    '''Abstract class for the two nucleotide coalescence system.

    Implements the basic state space exploration functionality, but
    leaves it to sub-classes to specify the actual system.

    The transitions in the system are specified in self.transitions
    that is a list of lists of transitions, where self.transitions[i]
    contains the transitions that take i+1 tokens as input.  The
    transitions consists of pairs, where the first element is a
    transition label (like "R" for recombination) and the second is
    the function that translates the pre-set into the post-set.

    For convenience, the pre-set is give as just the tokens (so not a
    set but a normal list of argumennts).

    The post-set should always be a set or a list of sets, since it is
    otherwise hard to figure out if the transition produces single
    tokens, multiple tokens or no tokens.

    If a set is returned, it is interpreted as the post set of firing
    the transition.  If a list (of sets) is returned, it is
    interpreted as if the transition can fire and
    non-deterministically produce different post sets; the lists
    should contain all those posts sets.
    
    The one exeception is if the transition returns None.  This is
    interpreted as a guard violation and the transition update is
    aborted.

    If the states are provided in the constructor (by default they are
    not) this is interpreted as the fixed state space and only
    transitions will be computed.  If the state space exploration
    encounters a state not in the provided list, it is considered an
    error and the state space exploration will be aborted.
    '''

    def __init__(self, species):
        self.species = species
        self.transitions = []
        self.state_numbers = None
        self.compute_states = True

    def successors(self, state):
        '''Calculate all successors of "state".

        Generates all successors of "state" and lets you iterate over the
        edges from "state" to tis successors.  Each generated value is a
        pair of transition type and state, where transition type is either
        "C" for a coalescence event or "R" for a recombination event.
        '''

        L = list(state)

        ttype,tfunc = self.transitions[0][0]
        for token in L:
            pre = iset([token])
            tproduct = tfunc(token)
            if tproduct:
                post = tproduct
                new_state = state.difference(pre).union(post)
                yield ttype, new_state

        ttype,tfunc = self.transitions[1][0]
        for i in xrange(len(L)):
            t1 = L[i]
            for j in xrange(i):
                t2 = L[j]
                
                pre = iset([t1,t2])

                tproduct = tfunc(t1,t2)
                post = tproduct
                if post == None:
                    continue
                new_state = state.difference(pre).union(post)
                yield ttype, new_state


    def compute_state_space(self):
        '''Computes the CTMC system for "species".'''
        initial_states = self.initial_state()

        seen = set(initial_states)
        unprocessed = initial_states
        self.state_numbers = {}
        for i,s in enumerate(initial_states):
            self.state_numbers[s] = i
        edges = []

        while unprocessed:
            s = unprocessed.pop()
            n = self.state_numbers[s]
            for t,ss in self.successors(s):
                assert s != ss, "We don't like self-loops!"

                if ss not in self.state_numbers:
                    assert self.compute_states, "Unknown state encountered while we only expect known states!"
                    self.state_numbers[ss] = len(self.state_numbers)

                if ss not in seen:
                    unprocessed.append(ss)
                    seen.add(ss)

                m = self.state_numbers[ss]
                edges.append((n,t,m))

        return self.state_numbers, edges



class BasicCoalSystem(CoalSystem):
    '''The basic two-nucleotide coalescence system.'''

    def recombination(self, token):
        '''Compute the tokens we get from a recombination on "token".
        
        Returns None if the recombination would just return
        "token" again, so we avoid returning "empty" tokens.
        '''
        _, nucs = token
        left, right = nucs
        if not (left and right): return None
        return iset([(0,(left,iset())), (0, (iset(),right))])

    def coalesce(self, token1, token2):
        '''Construct a new token by coalescening "token1" and "token2".'''
        _, nuc1 = token1
        _, nuc2 = token2
        left1, right1 = nuc1
        left2, right2 = nuc2
        left, right = left1.union(left2), right1.union(right2)
        return iset([(0, (left, right))])

    def initial_state(self):
        '''Build the initial state for this system.

        This doesn't necessarily mean that there is only a single
        initial state, but that doesn't matter much since we just need
        a state in an initial connected component for this to work...
        '''
        return [iset([(0,(iset([s]),iset([s]))) for s in self.species])]

    def __init__(self, species):
        CoalSystem.__init__(self, species)
        self.transitions = [[('R',self.recombination)], [('C',self.coalesce)]]

class SeperatedPopulationCoalSystem(CoalSystem):
    def recombination(self, token):
        '''Compute the tokens we get from a recombination on "token".
        
        Returns None if the recombination would just return
        "token" again, so we avoid returning "empty" tokens.
        '''
        pop, nucs = token
        left, right = nucs
        if not (left and right): return None # abort transition...
        return iset([(pop,(left,iset())),
                     (pop,(iset(),right))])

    def coalesce(self, token1, token2):
        '''Construct a new token by coalescening "token1" and "token2".'''
        pop1, nuc1 = token1
        pop2, nuc2 = token2
        if pop1 != pop2: return None # abort transition...

        left1, right1 = nuc1
        left2, right2 = nuc2
        left, right = left1.union(left2), right1.union(right2)
        return iset([(pop1,(left, right))])

    def initial_state(self):
        '''Build the initial state for this system.

        This doesn't necessarily mean that there is only a single
        initial state, but that doesn't matter much since we just need
        a state in an initial connected component for this to work...
        '''
        return self.init

    def __init__(self, species, initial_states=None):
        CoalSystem.__init__(self, species)
        self.transitions = [[('R',self.recombination)], [('C',self.coalesce)]]
        if initial_states == None:
            self.init = [iset([(s,(iset([s]),iset([s]))) for s in self.species])]
        else:
            self.init = initial_states

if __name__ == "__main__":
    pass

