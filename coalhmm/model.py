from scipy import identity, matrix, newaxis, zeros, array, int32, allclose, isnan
from scipy.linalg import expm
import scipy.weave as weave
iset = frozenset
from itertools import izip

from intervals import enumerate_all_transitions
from statespace_generator import SeperatedPopulationCoalSystem
from scc import SCCGraph, EpochSeperatedSCCGraph
from tree import make_tree
from emission_matrix import build_emission_matrix

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
    def __init__(self, nleaves, nbreakpoints, G, mappings=[]):
        self.nleaves = nleaves
        self.nbreakpoints = nbreakpoints
        self.G = G
        epoch_sizes = self.G.getEpochSizes()
        self.all_sizes = []
        for e in xrange(len(nbreakpoints)):
            self.all_sizes.extend([epoch_sizes[e]] * nbreakpoints[e])
        # Build a list of all paths through the SCC graph
        paths = []
        for S in G.all_paths():
            paths.append(S[:])

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
            assert 0 < b <= len(self.components_flat)
            self.components_flat[a:b] = array(c, dtype=int32)
            component_starts.append(a)
            component_ends.append(b)
        assert all(component_ends[i] == component_starts[i+1] for i in range(len(component_ends)-1))

        # Build all distributions of the paths over our intervals
        paths_final = []
        tree_map = {}
        paths_indices = []
        npaths = 0
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
            paths_final.extend(path_as_offsets)
            npaths += 1

            ta = make_tree(G, s, 0)
            tb = make_tree(G, s, 1)
            a = tree_map.setdefault(ta, len(tree_map))
            b = tree_map.setdefault(tb, len(tree_map))
            paths_indices.append(a)
            paths_indices.append(b)

        self.tree_map = tree_map
        self.ntrees = len(tree_map)
        self.paths_final_indices = array(paths_indices, dtype=int32)
        self.paths_final = array(paths_final, dtype=int32)
        self.npaths = npaths
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
            Q = zeros((n_states, n_states))
            for (src, transition, dst) in transitions:
                Q[src,dst] = rates[transition]
            for i in xrange(n_states):
                row = Q[i, :]
                Q[i,i] = -row.sum()
            return matrix(Q)

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
            assert (P > 0).any()
            assert not (isnan(P).any())
            e = in_epoch[j]
            if in_epoch[j+1] != e:
                fromSize = all_sizes[j]
                toSize = all_sizes[j+1]
                X = self.projMatrix(fromSize, toSize, mappings[e])
                P = P * X
            assert P.shape == (all_sizes[j], all_sizes[j+1])
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
        P_i_offsets = P_i_sizes.cumsum(dtype=int32)
        Pstart = zeros(P_i_sizes.sum())
        for i, P_i in enumerate(Ps):
            a = P_i_offsets[i]
            b = P_i_offsets[i+1]
            Pstart[a:b] = P_i[:]
            assert P_i.shape[0] == b-a
            assert 0 < b <= len(Pstart)
        pi_offsets = array([1]+all_sizes, dtype=int32).cumsum(dtype=int32)
        pi_buffer = zeros(sum(all_sizes)+1)
        components_flat = self.components_flat

        joint_prob_code = """
        #define myassert(c) do { if (!(c)) _myassert(#c, __LINE__, __FILE__); } while(0)
        void _myassert(const char *cond_str, int line, const char *filename)
        {
            printf("'%s' failed @ %s:%i\\n", cond_str, filename, line);
            exit(1);
        }
        double joint_prob(
            int nintervals,
            double *pi_buffer, int *pi_offsets,
            double *Pstart, int *P_i_offsets,
            int *components_flat, int *sizes,
            int *path_as_offsets)
        {
            int Sp = 0, Sc = 0;
            double *pi_curr = 0;
            double *pi_prev = pi_buffer;
            myassert(path_as_offsets[0] == 0);
            myassert(path_as_offsets[1] == 1);
            pi_prev[0] = 1.0;
            for (int i = 0; i < nintervals-1; i++)
            {
                int pa = path_as_offsets[2*i     + 0];
                int pb = path_as_offsets[2*i     + 1];
                int ca = path_as_offsets[2*(i+1) + 0];
                int cb = path_as_offsets[2*(i+1) + 1];
                double *P_i = Pstart + P_i_offsets[i];
                pi_curr = pi_buffer + pi_offsets[i+1];
                Sp = sizes[i]; Sc = sizes[i+1];
                for (int t = pa; t < pb; t++)
                {
                    int ct = components_flat[t];
                    // myassert(0 <= ct && ct < Sp);
                    for (int s = ca; s < cb; s++)
                    {
                        int cs = components_flat[s];
                        // myassert(0 <= cs && cs < Sc);
                        pi_curr[cs] += P_i[ct*Sc + cs] * pi_prev[ct];
                    }
                }
                pi_prev = pi_curr;
            }
            double res = 0.0;
            for (int i = 0; i < Sc; i++)
                res += pi_curr[i];
            return res;
        }
        """

        call_join_prob_code = """
        double total_joint = 0.0;
        for (int i = 0; i < npaths; i++)
        {
            int *path_as_offsets = paths_final + (i*2*nintervals);
            for (int j = 0; j < len_pi_buffer; j++)
                pi_buffer[j] = 0;
            double joint = joint_prob(
                    nintervals,
                    pi_buffer, pi_offsets,
                    Pstart, P_i_offsets,
                    components_flat, all_sizes,
                    path_as_offsets);
            int a = paths_final_indices[i*2 + 0];
            int b = paths_final_indices[i*2 + 1];
            J2(a, b) += joint;
            total_joint += joint;
        }
        return_val = total_joint;
        """
        all_sizes = array(all_sizes, dtype=int32)
        ntrees = len(tmap)
        nintervals = sum(self.nbreakpoints)
        npaths = self.npaths
        paths_final = self.paths_final
        assert len(paths_final) == npaths*2*nintervals
        paths_final_indices = self.paths_final_indices
        len_pi_buffer = len(pi_buffer)
        J = zeros((ntrees,ntrees))
        total_joint = weave.inline(call_join_prob_code,
                arg_names=['nintervals', 'npaths',
                    'paths_final', 'paths_final_indices',
                    'pi_buffer', 'pi_offsets', 'len_pi_buffer',
                    'components_flat', 'Pstart', 'P_i_offsets',
                    'J', 'all_sizes'],
                support_code=joint_prob_code,
                compiler='gcc')

        assert abs(total_joint - 1.0) < 0.0001
        # The starting probabilities are equal to the row-sums of J
        pi = J.sum(axis=1)
        assert abs(pi.sum() - 1.0) < 0.0001
        # The transitions have to be normalized
        T = J/pi[:,newaxis]
        assert all(abs(row.sum() - 1.0) < 0.0001 for row in T)

        Em = build_emission_matrix(tmap.keys(), tmap, col_map,\
                self.nleaves, breakpoints, in_epoch, theta, Qs, G, all_rates)
        # If we have a N column, each row will sum to 2 rather than 1
        if allclose(Em[:,-1].sum(), len(Em[:,-1])):
            assert all(allclose(row.sum(), 2.0) for row in Em)
        else:
            assert all(allclose(row.sum(), 1.0) for row in Em)
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
