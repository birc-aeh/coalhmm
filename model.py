from scipy import *
from scipy.linalg import expm
from sets import ImmutableSet as iset

from intervals import *
from statespace_generator import BasicCoalSystem
from scc import SCCGraph, EpochSeperatedSCCGraph
from tree import *
from emission_matrix import *
from time_plot import *

def prettify_state(s):
    """Convert a coal system state to something nicer.
    
    example:
      iset([(iset([3]), iset([3])),
            (iset([1]), iset([1])),
            (iset([2]), iset([2]))])
    to:
      {3, 1, 2}, {3, 1, 2}
    """
    def f(s, side, d):
        if d == 0:
            tmp = [f(sub, side, d+1) for sub in s if len(s) > 0]
            return ", ".join(sorted([x for x in tmp if x.strip() != ""]))
        elif d == 1:
            return "".join(sorted([str(x) for x in s[side] if len(s[side]) > 0]))
    return "{" + f(s, 0, 0) + "}, {" + f(s, 1, 0) + "}"


class Model:
    def __init__(self, nleaves, nintervals, G=None):
        self.nleaves = nleaves
        self.nintervals = nintervals
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
        for s in enumerate_all_transitions(paths, nintervals):
            paths_final.append(s)
            ta = make_tree(G, s, 0)
            if ta not in tree_map:
                tree_map[ta] = len(tree_map)
        self.tree_map = tree_map
        self.ntrees = len(tree_map)
        self.paths_final = paths_final

    def run(self, R, C, epoch_bps=None, mappings=[]):
        """Generates the parts needed for the HMM.
        Inititial state probabilities,
        Transition matrix, and
        Emmision matrix.
        Additionally returns the rate matrix used.
        """
        theta = 1 / C
        if epoch_bps == None:
            # TODO: choose better breakpoints?
            epoch_bps = [[0.0] + [c * theta for c in [.5,1,2,3,4]]]

        epoch_sizes = self.G.getEpochSizes()
        nepochs = len(epoch_sizes)
        assert sum(map(len,epoch_bps)) == self.nintervals + 1, \
                "We need n+1 breakpoints, for n intervals"
        assert len(mappings) == nepochs - 1, \
                "We need n-1 mappings for n epochs"

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
            Q = genRateMatrix(len(V), E, C=C, R=R)
            epochQ.append(Q)
            nbps = len(epoch_bps[e])
            Qs = Qs + [Q] * nbps
            in_epoch = in_epoch + [e] * nbps

        projections = []
        for j in xrange(len(breakpoints)-1):
            dt = breakpoints[j+1] - breakpoints[j]
            P.append(expm(Qs[j]*dt))
            if j > 1 and in_epoch[j-1] != in_epoch[j]:
                V, E = G.originalGraph(in_epoch[j-1])
                proj = arange(len(V))
                for a, b in mappings[in_epoch[j-1]]:
                    proj[a] = b
                print "mapping", j, proj
                print "from size", epoch_sizes[in_epoch[j-1]], "to", epoch_sizes[in_epoch[j]]
            else:
                V, E = G.originalGraph(in_epoch[j])
                proj = arange(len(V))
            projections.append(proj)
        print projections
        print [x.shape for x in P]
        print [x.shape for x in Qs]

        # Calculate the joint probability for a path through the graph
        # (The path must have an entry for each time interval)
        def joint_prob(sizes, G, path):
            # An extra state is added, as all paths should start in state 0,
            # meaning all species are seperate.
            all_seperate = 0
            component_path = [[all_seperate]] + [G.all_states(e, p) for (e,p) in path]
            print path
            lenV = sizes[0]
            pi_prev = zeros(lenV)
            pi_prev[all_seperate] = 1.0
            for i in xrange(len(component_path)):
                epoch = path[i][0]
                lenV = sizes[epoch]
                P_i = P[i]
                #TODO: matrices are not aligned?
                #pi_prev = projections[i] * pi_prev
                proj = projections[i]
                pi_curr = zeros(lenV)
                #print pi_prev, sum(pi_prev)
                #print i, in_epoch[i], proj
                for s in component_path[i+1]:
                    for x in component_path[i]:
                        print "asdasd", i, s, x, P_i.shape, P_i[proj[x], s] * pi_prev[x]
                        pi_curr[s] += P_i[x, s] * pi_prev[x]
                pi_prev = pi_curr
            return sum(pi_curr)

        ntrees = len(tmap)
        J = zeros((ntrees,ntrees))
        total_joint = 0.0
        for p in self.paths_final:
            joint = joint_prob(epoch_sizes, G, p)
            total_joint += joint
            a = tmap[make_tree(G, p, 0)]
            b = tmap[make_tree(G, p, 1)]
            J[a, b] = joint
        # TODO: reasonable epsilon?
        assert abs(total_joint - 1.0) < 0.0001

        # The starting probabilities is equal to the row-sums of J
        pi = sum(J, axis=0)
        # The transitions have to be normalized
        T = J/pi
        return pi, T, Em, Q[0] #TODO: don't return Q?

    def write_dot(self, filename):
        f = file(filename, "w")
        f.write("digraph {\n")
        for a in xrange(len(self.G.E)):
            f.write('%i [label="%s"]\n' % (a, prettify_state(self.G.state(a))))
            for b in self.G.E[a]:
                f.write("%i -> %i\n" % (a, b))
        f.write("}\n")
        f.close()

