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


x = BasicCoalSystem([0,1,2])
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
#for t in unique_topologies:
#    print tree_to_newick(t)


I = 5
interval_times = [0.0, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1]
def jukes_cantor(dt):
    gamma = 0.25 + 0.75 * exp(-4*dt)
    delta = 0.25 - 0.25 * exp(-4*dt)
    return matrix( [[gamma, delta, delta, delta],
                    [delta, gamma, delta, delta],
                    [delta, delta, gamma, delta],
                    [delta, delta, delta, gamma]])
  

S = [[None] * I for i in xrange(I)]
for i in xrange(I):
    for j in xrange(i, I):
        dt = interval_times[j]-interval_times[i]
        S[i][j] = jukes_cantor(dt)

for i in xrange(I):
    for j in xrange(i+1, I):
        #print S[i][j]
        pass

theta = 20000.0 * 20 * 1e-9

def D(i, j):
    t = interval_times
    def m(i):
        if i + 1 == len(t):
            dt = 0.0
            a = 0.0
        else:
            dt = t[i+1] - t[i]
            a = exp(-dt/theta)
        return theta - (dt * a)/(1 - a)
    dt_i = t[i+1] - t[i]
    return t[j] - t[i+1] + m(j) + dt_i - m(i)


test_tree = (4, iset([iset([1]), (2, iset([iset([0]), iset([2])]))]))

def first(s):
    for x in s:
        return x

to_char = ['A', 'C', 'G', 'T']
def emission_test(tree, S, cols):
    def visit(t):
        if isinstance(t, iset): # leaf
            res = array([0.0, 0.0, 0.0, 0.0])
            res[cols[first(t)]] = 1.0
            return 0, res
        else: # node
            j, children = t
            prob_j = ones(4)
            for arr in map(visit, children):
                i, prob_i = arr
                D_ij = D(i, j)
                S_ij = S[i][j]
                prob = zeros(4)
                for y in xrange(4):
                    for x in xrange(4):
                        prob[y] += D_ij * prob_i[x] * S_ij[x, y] 
                prob_j = prob_j * prob
            return j, prob_j
    return array([0.25, 0.25, 0.25, 0.25]) * visit(tree)[1]

print tree_to_newick(test_tree)
print emission_test(test_tree, S, [0,0,1])

