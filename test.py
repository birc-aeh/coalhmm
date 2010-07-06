from scipy import *
from scipy.linalg import expm
from sets import ImmutableSet as iset

from intervals import *
from statespace_generator import BasicCoalSystem
from scc import build_scc, SCCGraph
from tree import *
from emission_matrix import *

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
            return ", ".join([x for x in tmp if x.strip() != ""])
        elif d == 1:
            return "".join(sorted([str(x) for x in s[side] if len(s[side]) > 0]))
    return "{" + f(s, 0, 0) + "}, {" + f(s, 1, 0) + "}"


class Model:
    def __init__(self):
        pass

    def prepare(self, n):
        # Build transition system
        x = BasicCoalSystem(range(n))
        states, edges = x.compute_state_space()
        # Build SCC graph and do transitive closure
        G = SCCGraph(states, edges)
        G.add_transitive_edges()
        self.G = G
        # Build a list of all paths through the SCC graph
        paths = []
        def dfs(a, S, E):
            S.append(a)
            if len(E[a]) == 0:
                paths.append(S[:])
            for b in E[a]:
                dfs(b, S, E)
            S.pop()
        dfs(len(G.V)-1, [], G.E)
        # Build all distributions of the paths over our intervals
        paths_final = []
        tree_map = {}
        for s in enumerate_all_transitions(paths, 5):
            paths_final.append(s)
            ta = make_tree(G, s, 0)
            if ta not in tree_map:
                tree_map[ta] = len(tree_map)
            tb = make_tree(G, s, 1)
        self.tree_map = tree_map
        self.ntrees = len(tree_map)
        self.paths_final = paths_final

        self.ready = True

    def run(self, R, C, interval_times=None):
        if not self.ready:
            return

        theta = 1 / C
        if interval_times == None:
            interval_times = [0.0] + [c * theta for c in [.5,1,2,3,4]]

        tmap = self.tree_map
        G = self.G
        Em = build_emission_matrix(tmap.keys(), tmap, 3, interval_times, theta)

        def genRateMatrix(states,edges,**mapping):
            def f(t):
                return mapping[t]

            n_states = len(states)
            M = zeros((n_states, n_states))
            M.shape = n_states, n_states
            for (a,t,b) in edges:
                M[a,b] = f(t)
            for i in xrange(n_states):
                row = M[i, :]
                M[i,i] = -sum(row)

            M = matrix(M)
            return M

        Ps = []
        V, E = G.originalGraph()
        Q = genRateMatrix(V, E, C=C, R=R)
        for i in xrange(len(interval_times)-1):
            dt = interval_times[i+1] - interval_times[i]
            Ps.append(expm(Q*dt))

        def joint_prob(V, G, path):
            component_path = [[0]] + [G.all_states(p) for p in path]
            pi_prev = zeros(len(V))
            pi_prev[0] = 1.0
            for i in xrange(len(path)):
                pi_curr = zeros(len(V))
                P_i = Ps[i]
                for s in component_path[i+1]:
                    for x in component_path[i]:
                        pi_curr[s] += pi_prev[x] * P_i[x,s]
                pi_prev = pi_curr
            return sum(pi_curr)

        ntrees = len(tmap)
        J = zeros((ntrees,ntrees))
        total_joint = 0.0
        for p in self.paths_final:
            joint = joint_prob(V, G, p)
            #print p,">",joint
            total_joint += joint
            ta = tmap[make_tree(G, p, 0)]
            tb = tmap[make_tree(G, p, 1)]
            J[ta, tb] = joint
        # TODO: reasonable epsilon?
        assert abs(total_joint - 1.0) < 0.0001

        pi = sum(J, axis=0)
        T = J/pi
        return pi, T, Em

    def write_dot(self, filename):
        if not self.ready:
            return
        f = file(filename, "w")
        f.write("digraph {\n")
        for a in xrange(len(self.G.E)):
            f.write('%i [label="%s"]\n' % (a, prettify_state(self.G.state(a))))
            for b in self.G.E[a]:
                f.write("%i -> %i\n" % (a, b))
        f.write("}\n")
        f.close()


theta = 20000.0 * 20 * 1e-9
rho = 20000.0 * 0.01 / 1e6

model = Model()
model.prepare(2)

E, T, pi = model.run(rho * (1.0 / theta), 1.0 / theta)
#print "lala"
#E, T, pi = model.run(rho * (2.0 / theta), 2.0 / theta)

