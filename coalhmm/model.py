from scipy import identity, asmatrix, newaxis, zeros, array, int32
from scipy.linalg import expm
import scipy.weave as weave
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
        epoch_sizes = self.G.getEpochSizes()
        self.all_sizes = []
        for e in xrange(len(nbreakpoints)):
            self.all_sizes.extend([epoch_sizes[e]] * nbreakpoints[e])
        # Build a list of all paths through the SCC graph
        paths = []
        for S in G.all_paths():
            paths.append(S[:])
        #print "All paths"
        #for p in paths:
        #    print ["e%i_c%i" % (e,c) for e,c in p]

        epoch_sizes = G.getEpochSizes()
        component_index = dict()
        components = [(0,)]
        for e,esize in enumerate(epoch_sizes):
            for c in xrange(len(G.G[e].V)):
                component = array(G.all_states(e,c))
                component_idx = len(components)
                components.append(component)
                component_index[(e,c)] = component_idx

        self.components_flat = zeros(sum(self.all_sizes)+1, dtype=int32)
        component_starts = []
        component_ends = []
        offset = 0
        for c in components:
            a = offset
            b = offset + len(c)
            offset = b
            assert self.components_flat[a:b].sum() == 0
            self.components_flat[a:b] = array(c, dtype=int32)
            component_starts.append(a)
            component_ends.append(b)

        # Build all distributions of the paths over our intervals
        paths_final = []
        tree_map = {}
        paths_indices = []
        for s in enumerate_all_transitions(paths, nbreakpoints):
            # FIXME: instead of removing the first component in the path,
            # we shouldn't have it there to begin with...
            s = s[1:]
            cpath = (0,)+tuple(component_index[(e,p)] for e,p in s)
            path_as_offsets = []
            for ci in cpath:
                path_as_offsets.append(component_starts[ci])
                path_as_offsets.append(component_ends[ci])
            path_as_offsets = array(path_as_offsets, dtype=int32)
            paths_final.append(path_as_offsets)

            ta = make_tree(G, s, 0)
            tb = make_tree(G, s, 1)
            a = tree_map.setdefault(ta, len(tree_map))
            b = tree_map.setdefault(tb, len(tree_map))
            paths_indices.append((a, b))

        self.tree_map = tree_map
        self.ntrees = len(tree_map)
        self.paths_final_indices = paths_indices
        self.paths_final = paths_final
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
        assert len(Qs) == len(breakpoints)

        all_sizes = self.all_sizes
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
            Ps.append(array(P).flatten())
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

        P_i_sizes = array([0]+[P_i.shape[0] for P_i in Ps], dtype=int32)
        self.P_i_offsets = P_i_sizes.cumsum()
        self.Pstart = zeros(P_i_sizes.sum())
        for i, P_i in enumerate(Ps):
            a = self.P_i_offsets[i]
            b = self.P_i_offsets[i+1]
            self.Pstart[a:b] = P_i[:]
        self.pi_offsets = array([1]+all_sizes, dtype=int32).cumsum()
        pi_buffer = zeros(sum(all_sizes)+1)
        def joint_prob_cached(sizes, path_as_offsets):
            code = """
            int Sp = 0, Sc = 0;
            double *pi_curr = 0;
            double *pi_prev = pi_buffer;
            pi_prev[0] = 1.0;
            for (int i = 0; i < PATH_LEN-1; i++)
            {
                int pa = path_as_offsets[2*i     + 0];
                int pb = path_as_offsets[2*i     + 1];
                int ca = path_as_offsets[2*(i+1) + 0];
                int cb = path_as_offsets[2*(i+1) + 1];
                double *P_i = Pstart + P_i_offsets[i];
                pi_curr = pi_buffer + pi_offsets[i];
                Sp = sizes[i]; Sc = sizes[i+1];
                for (int t = pa; t < pb; t++)
                {
                    int ct = components_flat[t];
                    for (int s = ca; s < cb; s++)
                    {
                        int cs = components_flat[s];
                        pi_curr[cs] += P_i[ct*Sc + cs] * pi_prev[ct];
                    }
                }
                pi_prev = pi_curr;
            }
            double res = 0.0;
            for (int i = 0; i < Sc; i++)
                res += pi_curr[i];
            return_val = res;
            """
            components_flat = self.components_flat
            pi_buffer[:] = 0.0
            pi_offsets = self.pi_offsets
            Pstart = self.Pstart
            P_i_offsets = self.P_i_offsets
            PATH_LEN = len(path_as_offsets)/2

            return weave.inline(code,
                    ['sizes', 'pi_buffer', 'pi_offsets',
                        'Pstart', 'P_i_offsets',
                        'PATH_LEN',
                        'components_flat', 'path_as_offsets'],
                    compiler='gcc')

        all_sizes = array(all_sizes, dtype=int32)
        ntrees = len(tmap)
        J = zeros((ntrees,ntrees))
        total_joint = 0.0
        for p, (a,b) in izip(
                self.paths_final,
                self.paths_final_indices):
            joint = joint_prob_cached(all_sizes, p)
            assert joint >= 0
            total_joint += joint
            J[a, b] += joint

        # TODO: reasonable epsilon?
        assert abs(total_joint - 1.0) < 0.0001
        # The starting probabilities are equal to the row-sums of J
        pi = J.sum(axis=1)
        assert abs(pi.sum() - 1.0) < 0.0001
        # The transitions have to be normalized
        T = J/pi[:,newaxis]
        assert all(abs(row.sum() - 1.0) < 0.0001 for row in T)

        Em = build_emission_matrix(tmap.keys(), tmap, col_map,\
                self.nleaves, breakpoints, in_epoch, theta, Qs, G, all_rates)
        return pi, T, Em

def build_epoch_separated_scc(nleaves, mappings, migration=None, init=None):
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
    statespace = SeperatedPopulationCoalSystem(range(nleaves), initial_states=init, legal_migrations=migration[0])
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

def build_epoch_separated_model(nleaves, mappings, epoch_nbps, migration=None, init=None):
    '''Creates a model separated in to epochs.

>>> m = build_epoch_separated_model(2, [], [3])
>>> m = build_epoch_separated_model(2, [[0,0]], [2,3])'''
    assert len(mappings) == len(epoch_nbps) - 1
    mappings, G = build_epoch_seperated_scc(nleaves, mappings, migration, init)
    return Model(nleaves, epoch_nbps, G, mappings)
build_epoch_seperated_model = build_epoch_separated_model

if __name__ == "__main__":
    import doctest
    doctest.testmod()
