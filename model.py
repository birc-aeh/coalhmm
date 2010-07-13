from scipy import *
from scipy.linalg import expm
from sets import ImmutableSet as iset

from intervals import *
from statespace_generator import BasicCoalSystem, SeperatedPopulationCoalSystem
from scc import SCCGraph, EpochSeperatedSCCGraph
from tree import *
from emission_matrix import *

def prettify_state(s):
    """Convert a coal system state to something nicer.
    
    example:
>>> from sets import ImmutableSet as iset
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
        def paths_filler(S):
            paths.append(S[:])
        G.all_paths(paths_filler)
        # Build all distributions of the paths over our intervals
        paths_final = []
        tree_map = {}
        # We assume one less breakpoint in the first epoch, because we later
        #  have to add a state where everything is seperated.
        nbreakpoints[0] = nbreakpoints[0] - 1
        for s in enumerate_all_transitions(paths, nbreakpoints):
            paths_final.append(s)
            ta = make_tree(G, s, 0)
            tb = make_tree(G, s, 1)
            if ta not in tree_map:
                tree_map[ta] = len(tree_map)
            if tb not in tree_map:
                tree_map[tb] = len(tree_map)
        nbreakpoints[0] = nbreakpoints[0] + 1 # bump it back up again
        self.tree_map = tree_map
        self.ntrees = len(tree_map)
        self.paths_final = paths_final
        self.mappings = mappings

    def run(self, R, C, epoch_bps):
        """Generates the parts needed for the HMM.
        Inititial state probabilities,
        Transition matrix, and
        Emmision matrix.
        Additionally returns the rate matrix used.
        """
        epoch_sizes = self.G.getEpochSizes()
        nepochs = len(epoch_sizes)
        mappings = self.mappings
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
        # TODO: transition can use R/C in different epochs, emission can't
        theta = 1 / (sum(C)/len(C))

        tmap = self.tree_map
        G = self.G
        breakpoints = []
        for x in epoch_bps:
            for t in x:
                breakpoints.append(t)
        Em = build_emission_matrix(tmap.keys(), tmap, self.nleaves, breakpoints, theta)

        def genRateMatrix(n_states,edges,**mapping):
            def f(t):
                return mapping[t]
            M = asmatrix(zeros((n_states, n_states)))
            for (a,t,b) in edges:
                M[a,b] = f(t)
            for i in xrange(n_states):
                row = M[i, :]
                M[i,i] = -sum(row)
            return M

        P = []
        epochQ = []
        Qs = []
        in_epoch = []
        for e in xrange(len(epoch_bps)):
            V, E = G.originalGraph(e)
            Q = genRateMatrix(len(V), E, C=C[e], R=R[e])
            epochQ.append(Q)
            nbps = len(epoch_bps[e])
            Qs = Qs + [Q] * (nbps)
            in_epoch = in_epoch + [e] * nbps

        projections = []
        for j in xrange(len(breakpoints)-1):
            dt = breakpoints[j+1] - breakpoints[j]
            P.append(expm(Qs[j+1]*dt))
            if in_epoch[j+1] != in_epoch[j]:
                V, E = G.originalGraph(in_epoch[j])
                proj = arange(len(V))
                for a, b in mappings[in_epoch[j]]:
                    proj[a] = b
            else:
                V, E = G.originalGraph(in_epoch[j])
                proj = arange(len(V))
            projections.append(proj)
        assert len(P) == len(breakpoints) - 1
        assert len(Qs) == len(breakpoints)

        # Calculate the joint probability for a path through the graph
        # (The path must have an entry for each time interval)
        def joint_prob(sizes, G, path):
            # An extra state is added, as all paths should start in state 0,
            # meaning all species are seperate.
            all_seperate = 0
            component_path = [[all_seperate]] + [G.all_states(e, p) for (e,p) in path]
            in_epoch = [0] + [e for (e,p) in path]
            lenV = sizes[0]
            pi_prev = zeros(lenV)
            pi_prev[all_seperate] = 1.0
            for i in xrange(len(component_path)-1):
                P_i = P[i]
                proj = projections[i]
                pi_curr = zeros(sizes[in_epoch[i+1]])
                for s in component_path[i+1]:
                    for x in component_path[i]:
                        pi_curr[s] += P_i[proj[x], s] * pi_prev[x]
                pi_prev = pi_curr
            return sum(pi_curr)

        ntrees = len(tmap)
        J = zeros((ntrees,ntrees))
        total_joint = 0.0
        for p in self.paths_final:
            joint = joint_prob(epoch_sizes, G, p)
            total_joint += joint
            t1 = make_tree(G, p, 0)
            t2 = make_tree(G, p, 1)
            a = tmap[t1]
            b = tmap[t2]
            J[a, b] += joint

        # TODO: reasonable epsilon?
        assert abs(total_joint - 1.0) < 0.0001
        # The starting probabilities are equal to the row-sums of J
        pi = sum(J, axis=0)
        # The transitions have to be normalized
        T = J/pi
        return pi, T, Em

def build_epoch_seperated_scc(nleaves, mappings):
    '''Takes in a number of species, and a series of projections to build a
    SCC graph, seperated in to epochs.
    The projections (or mappings) consists of arrays - each species is given
    a number from 0 to nleaves-1, and looks in a projection at that index to
    find the new index.
    Merging species 0 and 1 out of three would be [0,0,1].
    
>>> proj, g = build_epoch_seperated_scc(3, [[0, 0, 0]])
>>> g.getEpochSizes()
[8, 203]
    '''
    statespace = SeperatedPopulationCoalSystem(range(nleaves))
    states, edges = statespace.compute_state_space()
    SCCs = [SCCGraph(states, edges, 0)]
    epoch = 1
    for mapping in mappings:
        init = [iset([(mapping[p], tok) for (p, tok) in x]) for x in states]
        statespace = SeperatedPopulationCoalSystem(range(nleaves), init)
        states, edges = statespace.compute_state_space()
        G = SCCGraph(states, edges, epoch)
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
                m = G.add_component_edge(i, old_state, i+1,new_state)
                tmp_proj.append(m)
        projections.append(tmp_proj)
    return projections, G

def build_simple_model(nleaves, bps):
    '''Creates a model using the simple statespace, where all leaves are
    considered as part of the same population.
    
>>> m = build_simple_model(2, 3)'''
    return Model(nleaves, [bps])

def build_epoch_seperated_model(nleaves, mappings, epoch_nbps):
    '''Creates a model seperated in to epochs.
    
>>> m = build_epoch_seperated_model(2, [], [3])
>>> m = build_epoch_seperated_model(2, [[0,0]], [2,3])'''
    assert len(mappings) == len(epoch_nbps) - 1
    mappings, G = build_epoch_seperated_scc(nleaves, mappings)
    return Model(nleaves, epoch_nbps, G, mappings)

if __name__ == "__main__":
    import doctest
    doctest.testmod()
