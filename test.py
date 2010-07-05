from statespace_generator import BasicCoalSystem
from scc import build_scc
from scipy import *
from scipy.linalg import expm
from sets import ImmutableSet as iset
from intervals import *

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

def prettify_set(s):
    def f(s, side, d):
        if d == 0:
            tmp = [f(sub, side, d+1) for sub in s if len(s) > 0]
            return ", ".join([x for x in tmp if x.strip() != ""])
        elif d == 1:
            return "".join(sorted([str(x) for x in s[side] if len(s[side]) > 0]))
    return "{" + f(s, 0, 0) + "}, {" + f(s, 1, 0) + "}"


x = BasicCoalSystem([1,2,3])
states, edges = x.compute_state_space()
states_rev = {}
for k, v in states.iteritems():
    states_rev[v] = k
#print sorted(states.values())
edges2 = [[] for i in xrange(len(states))]
for (a,t,b) in edges:
    edges2[a].append(b)
V, E = build_scc(states, edges2)
print '\n\n'.join(['\n'.join(map(prettify_set,[states_rev[x] for x in sorted(c)])) for c in V])
print "len:", len(V)
#print E
#print genRateMatrix(states,edges,C=1.0,R=1.0e-4)

#for n,t,m in edges:
#    print '%s -> %s -> %s' % (str(n),t,str(m))



for i in xrange(len(states_rev)):
    print '%d : %s' % (i,prettify_set(states_rev[i]))


def transitive_edges(a, S, newE):
    for v in S:
        newE[v].add(a)
    S.append(a)
    for b in E[a]:
        transitive_edges(b, S, newE)
    S.pop()
E2 = [set() for x in E]
transitive_edges(len(V)-1, [], E2)
print E
print E2

if True:
    f = file("graph", "w")
    f.write("digraph {\n")
    for a in xrange(len(E)):
        f.write('%i [label="%s"]\n' % (a, prettify_set(states_rev[V[a][0]])))
        for b in E2[a]:
            f.write("%i -> %i\n" % (a, b))
    f.write("}\n")
    f.close()

def dfs(a, S, E):
    S.append(a)
    if len(E[a]) == 0:
        print "path:", S
        print_all_distributions(S, 6)
    for b in E[a]:
        dfs(b, S, E)
    S.pop()

dfs(len(V)-1, [], E2)
