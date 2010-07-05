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
#print '\n\n'.join(['\n'.join(map(prettify_set,[states_rev[x] for x in sorted(c)])) for c in V])
#print "len:", len(V)
#print E
#print genRateMatrix(states,edges,C=1.0,R=1.0e-4)

#for n,t,m in edges:
#    print '%s -> %s -> %s' % (str(n),t,str(m))



#for i in xrange(len(states_rev)):
#    print '%d : %s' % (i,prettify_set(states_rev[i]))


def transitive_edges(a, S, newE):
    for v in S:
        newE[v].add(a)
    S.append(a)
    for b in E[a]:
        transitive_edges(b, S, newE)
    S.pop()
E2 = [set() for x in E]
transitive_edges(len(V)-1, [], E2)
#print E
#print E2

if True:
    f = file("graph", "w")
    f.write("digraph {\n")
    for a in xrange(len(E)):
        f.write('%i [label="%s"]\n' % (a, prettify_set(states_rev[V[a][0]])))
        for b in E2[a]:
            f.write("%i -> %i\n" % (a, b))
    f.write("}\n")
    f.close()


def project(s, side):
    return iset([x[side] for x in s if x[side] != iset()])

def make_tree(s):
    tree = None
    initial = project(states_rev[0], 0)
    used = set()
    for i in xrange(len(s)):
        if i == 0:
            B = project(states_rev[V[s[i]][0]], 0)
            if len(initial) != len(B):
                joined_from_the_start = initial - B
                tree = (0, [x for x in joined_from_the_start])
                for x in joined_from_the_start:
                    used.add(x)
        elif s[i] != s[i-1]:
            A = project(states_rev[V[s[i-1]][0]], 0)
            B = project(states_rev[V[s[i]][0]], 0)
            #print A
            #print B
            #print "merged:",(A-B)
            A = iset([x for x in A if len(x) == 1])
            joined = A-B
            if len(joined) == 0:
                continue
            for x in joined:
                used.add(x)
            if tree == None:
                tree = (i-1, [x for x in joined])
            else:
                tree = (i-1, [tree] + [x for x in joined])
    if len(initial) == len(used):
        return tree
    if len(used) == 0:
        rest_joined_at = 0
    else:
        rest_joined_at = len(s)
    rest = [x for x in initial if not x in used]
    return (rest_joined_at, (tree and [tree] or []) + rest)

#test = [24, 23, 23, 23, 23, 23, 20, 15]#[224, 217, 216, 208, 130, 0]

def tree_to_newick(t):
    if isinstance(t, iset):
        return "".join([str(x) for x in t])
    else:
        return "(" + ", ".join(map(tree_to_newick, t[1])) + "):" + str(t[0])

def test(s):
    print s, "-->", tree_to_newick(make_tree(s))

def dfs(a, S, E):
    S.append(a)
    if len(E[a]) == 0:
        print " -- starting on path:", S
        do_on_all_distributions(S, 5, test)
    for b in E[a]:
        dfs(b, S, E)
    S.pop()

dfs(len(V)-1, [], E2)

#for i in test:
#    print '%d : %s' % (i,prettify_set(states_rev[V[i][0]]))

