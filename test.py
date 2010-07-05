from scipy import *
from scipy.linalg import expm
from sets import ImmutableSet as iset

from intervals import *
from statespace_generator import BasicCoalSystem
from scc import build_scc, SCCGraph
from tree import *

#print genRateMatrix(states,edges,C=1.0,R=1.0e-4)
def genRateMatrix(states,edges,**mapping):
    def f(t):
        return mapping[t]

    n_states = len(states)
    print n_states
    M = zeros((n_states, n_states))
    M.shape = n_states, n_states
    for (a,t,b) in edges:
        M[a,b] = f(t)
    for i in xrange(n_states):
        row = M[i, :]
        M[i,i] = -sum(row)

    M = matrix(M)
    return M

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


x = BasicCoalSystem([1,2,3])
states, edges = x.compute_state_space()

G = SCCGraph(states, edges)
G.add_transitive_edges()

if False:
    f = file("graph", "w")
    f.write("digraph {\n")
    for a in xrange(len(G.E)):
        f.write('%i [label="%s"]\n' % (a, prettify_state(states_rev[V[a][0]])))
        for b in G.E[a]:
            f.write("%i -> %i\n" % (a, b))
    f.write("}\n")
    f.close()

def set_filler(S):
    def f(s):
        S.add(make_tree(G, s))
    return f

unique_topologies = set()
def dfs(a, S, E):
    S.append(a)
    if len(E[a]) == 0:
        do_on_all_distributions(S, 5, set_filler(unique_topologies))
    for b in E[a]:
        dfs(b, S, E)
    S.pop()

dfs(len(G.V)-1, [], G.E)
for t in unique_topologies:
    print tree_to_newick(t)

