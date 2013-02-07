from scipy import identity, asmatrix
from scipy.linalg import expm
iset = frozenset
from itertools import izip

from intervals import *
from statespace_generator import BasicCoalSystem, SeperatedPopulationCoalSystem, IM
from scc import SCCGraph, EpochSeperatedSCCGraph
from tree import *
from emission_matrix import *

def prettify_state(s):
    """Convert a coal system state to something nicer.
    
    example:
>>> iset = frozenset
>>> prettify_state(iset([(iset([3]), iset([3])),  (iset([1]), iset([1])),  (iset([2]), iset([2]))]))
'{1, 2, 3}, {1, 2, 3}'
    """
    def f(s, side, d):
        if d == 0:
            tmp = [f(sub, side, d+1) for sub in s if len(s) > 0]
            return ", ".join(sorted([x for x in tmp if x.strip() != ""]))
        elif d == 1:
            return "".join(sorted([str(x) for x in s[side] if len(s[side]) > 0]))
    return "{" + f(s, 0, 0) + "}, {" + f(s, 1, 0) + "}"


class Model:
    def __init__(self, nleaves, nbreakpoints, G=None, mappings=[]):
        self.nleaves = nleaves
        self.nbreakpoints = nbreakpoints
        if G == None:
            # Build transition system
            x = BasicCoalSystem(range(nleaves))
            states, edges = x.compute_state_space()
            # Build SCC graph and do transitive closure
            sG = SCCGraph(states, edges)
            sG.add_transitive_edges()
            G = EpochSeperatedSCCGraph()
            G.addSubGraph(sG)
        self.G = G
        # Build a list of all paths through the SCC graph
        paths = []
        for S in G.all_paths():
            paths.append(S[:])
        #print "All paths"
        #for p in paths:
        #    print ["e%i_c%i" % (e,c) for e,c in p]
        # Build all distributions of the paths over our intervals
        paths_final = []
        tree_map = {}
        paths_indices = []
        path_prefix_numbers = {}
        # For each path, p, paths_prefix_ids contains a list numbering the
        # prefixes of p. So the paths aaa and aab might be [0,1,2] and [0,1,3].
        # This is used later to do cache lookups without hasing a full path.
        paths_prefix_ids = []

        for s in enumerate_all_transitions(paths, nbreakpoints):
            # FIXME: instead of removing the first component in the path,
            # we shouldn't have it there to begin with...
            s = s[1:]
            paths_final.append(((0,),)+tuple(G.all_states(e,p) for e,p in s))
            ta = make_tree(G, s, 0)
            tb = make_tree(G, s, 1)
            a = tree_map.setdefault(ta, len(tree_map))
            b = tree_map.setdefault(tb, len(tree_map))
            paths_indices.append((a, b))
            ppn = path_prefix_numbers
            prefixes = [ppn.setdefault(s[:i+1], len(ppn)) for i in xrange(len(s))]
            paths_prefix_ids.append(prefixes)

        self.tree_map = tree_map
        self.ntrees = len(tree_map)
        # The paths are sorted by the prefix ids to give as much overlap as
        # possible between to consecutive paths (so we can utilize the cache
        # better)
        indices = sorted(range(len(paths_final)), key=lambda i: paths_prefix_ids[i])
        self.paths_final_indices = [paths_indices[i] for i in indices]
        self.paths_final = [paths_final[i] for i in indices]
        self.paths_prefix_ids = [array(paths_prefix_ids[i]) for i in indices]
        self.mappings = mappings

    def projMatrix(self, fromSize, toSize, mapping):
        res = zeros((fromSize, toSize))
        for a,b in mapping:
            res[a,b] = 1
        return matrix(res)

    def run(self, R, C, epoch_bps, M=None, col_map=None):
        """Generates the parts needed for the HMM.
        Inititial state probabilities,
        Transition matrix, and
        Emmision matrix.
        Additionally returns the rate matrix used.
        """
        epoch_sizes = self.G.getEpochSizes()
        nepochs = len(epoch_sizes)
        mappings = self.mappings
        nleaves = self.nleaves
        assert len(epoch_bps) == nepochs, \
                "We need breakpoints for every epoch"
        assert map(len,epoch_bps) == self.nbreakpoints, \
                "Wrong number of breakpoints"
        assert len(mappings) == nepochs - 1, \
                "We need n-1 mappings for n epochs"
        if not isinstance(R, list):
            assert not isinstance(C, list)
            R = [R] * nepochs
            C = [C] * nepochs
        assert isinstance(R, list) and isinstance(C, list)
        assert len(R) == len(C) == nepochs
        if not isinstance(R[0], list):
            assert not isinstance(C[0], list)
            R = [[r]*nleaves for r in R]
            C = [[c]*nleaves for c in C]
        else:
            assert map(len, R) == map(len, C)
        # TODO: transition can use R/C in different epochs, emission can't
        def avg(L):
            return sum(L)/len(L)
        theta = 1 / avg([avg(l) for l in C])

        # If no col_map is given, we assume that only A,C,G and T are possible
        if col_map == None:
            acgt = ['A', 'C', 'G', 'T']
            def i_to_cols(value, n):
                return map(lambda i: acgt[(value >> 2*i) & 3], xrange(n))
            col_map = dict()
            for i in xrange(1 << nleaves*2):
                col_map.setdefault(tuple(i_to_cols(i, nleaves)), len(col_map))

        if M == None:
            M = [identity(nleaves)] * nepochs

        tmap = self.tree_map
        G = self.G
        breakpoints = []
        for x in epoch_bps:
            for t in x:
                breakpoints.append(t)

        def genRateMatrix(n_states, transitions, rates):
            Q = asmatrix(zeros((n_states, n_states)))
            for (src, transition, dst) in transitions:
                Q[src,dst] = rates[transition]
            for i in xrange(n_states):
                row = Q[i, :]
                Q[i,i] = -row.sum()
            return Q

        Qs = []
        in_epoch = []
        all_sizes = []
        all_rates = []
        for e in xrange(len(epoch_bps)):
            V, E = G.originalGraph(e)
            rates = dict()
            pops = range(nleaves)
            for pa in pops:
                for pb in pops:
                    if pa == pb:
                        rates[('C',pa,pb)] = C[e][pa]
                        rates[('R',pa,pb)] = R[e][pa]
                    if pa != pb:
                        rates[('M',pa,pb)] = M[e][pa,pb]
            Q = genRateMatrix(len(V), E, rates)
            nbps = len(epoch_bps[e])
            Qs = Qs + [Q] * (nbps)
            all_rates = all_rates + [rates] * nbps
            in_epoch = in_epoch + [e] * nbps
            all_sizes = all_sizes + [epoch_sizes[e]] * nbps
        assert len(Qs) == len(breakpoints)

        Ps = []
        for j in xrange(len(breakpoints)-1):
            dt = breakpoints[j+1] - breakpoints[j]
            P = expm(Qs[j]*dt)
            assert dt >= 0
            assert (P >= 0).all()
            e = in_epoch[j]
            if in_epoch[j+1] != e:
                fromSize = all_sizes[j]
                toSize = all_sizes[j+1]
                X = self.projMatrix(fromSize, toSize, mappings[e])
                P = matrix(P) * matrix(X)
            Ps.append(P)
        assert len(Ps) == len(breakpoints) - 1

        # Calculate the joint probability for a path through the graph
        # (The path must have an entry for each time interval)
        # def joint_prob(sizes, path):
        #     # An extra state is added, as all paths should start in state 0,
        #     # meaning all species are seperate.
        #     pi_prev = zeros(sizes[0])
        #     pi_prev[0] = 1.0
        #     for i in xrange(len(path)-1):
        #         P_i = Ps[i]
        #         pi_curr = zeros(sizes[i+1])
        #         for s in path[i+1]:
        #             for x in path[i]:
        #                 pi_curr[s] += P_i[x, s] * pi_prev[x]
        #         pi_prev = pi_curr
        #     return sum(pi_curr)

        joint_prob_cache = {}
        def joint_prob_cached(sizes, path, path_prefix_ids):
            def jp(i):
                if i == 0:
                    res = zeros(sizes[0])
                    res[0] = 1.0
                    return res
                sub_path = path_prefix_ids[i-1]
                if sub_path in joint_prob_cache:
                    return joint_prob_cache[sub_path]
                P_i = Ps[i-1]
                pi_prev = jp(i-1)
                pi_curr = zeros(sizes[i])
                for s in path[i]:
                    for x in path[i-1]:
                        pi_curr[s] += P_i[x, s] * pi_prev[x]
                joint_prob_cache[sub_path] = pi_curr
                return pi_curr
            return sum(jp(len(path)-1))

        ntrees = len(tmap)
        J = zeros((ntrees,ntrees))
        total_joint = 0.0
        for p, pp, (a,b) in izip(
                self.paths_final,
                self.paths_prefix_ids,
                self.paths_final_indices):
            joint = joint_prob_cached(all_sizes, p, pp)
            assert joint >= 0
            total_joint += joint
            J[a, b] += joint

        # TODO: reasonable epsilon?
        assert abs(total_joint - 1.0) < 0.0001
        # The starting probabilities are equal to the row-sums of J
        pi = J.sum(axis=1)
        # The transitions have to be normalized
        T = J/pi

        Em = build_emission_matrix(tmap.keys(), tmap, col_map,\
                self.nleaves, breakpoints, in_epoch, theta, Qs, G, all_rates)
        return pi, T, Em

def build_epoch_separated_scc(nleaves, mappings, migration=None):
    '''Takes in a number of species, and a series of projections to build a
    SCC graph, seperated in to epochs.
    The projections (or mappings) consists of arrays - each species is given
    a number from 0 to nleaves-1, and looks in a projection at that index to
    find the new index.
    Merging species 0 and 1 out of three would be [0,0,1].
    
>>> proj, g = build_epoch_separated_scc(3, [[0, 0, 0]])
>>> g.getEpochSizes()
[8, 203]
    '''
    if not migration:
        migration = [None] * (1 + len(mappings))
    statespace = SeperatedPopulationCoalSystem(range(nleaves), legal_migrations=migration[0])
    states, edges = statespace.compute_state_space()
    SCCs = [SCCGraph(states, edges, 0, migration[0]!=None)]
    epoch = 1
    for mapping in mappings:
        init = [iset([(mapping[p], tok) for (p, tok) in x]) for x in states]
        statespace = SeperatedPopulationCoalSystem(range(nleaves), init, migration[epoch])
        states, edges = statespace.compute_state_space()
        G = SCCGraph(states, edges, epoch, migration[epoch]!=None)
        epoch += 1
        G.add_transitive_edges()
        SCCs.append(G)

    G = EpochSeperatedSCCGraph()
    G.addSubGraph(SCCs[0])
    projections = []
    for i, new in enumerate(SCCs[1:]):
        G.addSubGraph(new)
        old = SCCs[i] # note: i starts at 0, so 'new' is at i+1
        tmp_proj = []
        for c in xrange(len(old.V)):
            for v in old.all_states(c):
                old_state = old.states_rev[v]
                new_state = iset([(mappings[i][p], tok) for (p, tok) in old_state])
                m = G.add_component_edge(i, old_state, i+1, new_state)
                tmp_proj.append(m)
        projections.append(tmp_proj)
    return projections, G
build_epoch_seperated_scc = build_epoch_separated_scc

def build_simple_model(nleaves, bps):
    '''Creates a model using the simple statespace, where all leaves are
    considered as part of the same population.
    
>>> m = build_simple_model(2, 3)
>>> pi, T, E = m.run(0.1, 0.1, [[0.1, 0.2, 0.3]])
>>> T.shape, E.shape
((3, 3), (3, 16))'''
    return Model(nleaves, [bps])

def build_epoch_separated_model(nleaves, mappings, epoch_nbps, migration=None):
    '''Creates a model separated in to epochs.

>>> m = build_epoch_separated_model(2, [], [3])
>>> m = build_epoch_separated_model(2, [[0,0]], [2,3])'''
    assert len(mappings) == len(epoch_nbps) - 1
    mappings, G = build_epoch_seperated_scc(nleaves, mappings, migration)
    return Model(nleaves, epoch_nbps, G, mappings)
build_epoch_seperated_model = build_epoch_separated_model

if __name__ == "__main__":
    import doctest
    doctest.testmod()
