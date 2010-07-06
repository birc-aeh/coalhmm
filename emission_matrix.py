from scipy import *
from sets import ImmutableSet as iset

def first(s):
    """At the leaves of a tree, there is always isets of size one.
       This method is used to extract that member."""
    for x in s:
        return x

# a,b,c,d <=> a + (b << 2) + (c << 4) + (d << 6)
def cols_to_index(cols):
    i = 0
    s = 0
    for v in cols:
        s += v << i
        i += 2
    return s
def index_to_cols(value, n):
    return map(lambda i: (value >> 2*i) & 3, xrange(n))


def jukes_cantor(dt):
    gamma = 0.25 + 0.75 * exp(-4*dt)
    delta = 0.25 - 0.25 * exp(-4*dt)
    return matrix( [[gamma, delta, delta, delta],
                    [delta, gamma, delta, delta],
                    [delta, delta, gamma, delta],
                    [delta, delta, delta, gamma]])

def emission_test(tree, S, cols):
    def visit(t):
        if isinstance(t, iset): # leaf
            res = array([0.0, 0.0, 0.0, 0.0])
            res[cols[first(t)]] = 1.0
            return 0, res
        else: # node
            j, children = t
            prob_j = ones(4)
            for i, prob_i in map(visit, children):
                S_ij = S[i][j]
                prob = zeros(4)
                for y in xrange(4):
                    for x in xrange(4):
                        prob[y] += prob_i[x] * S_ij[x, y] 
                prob_j = prob_j * prob
            return j, prob_j
    return sum(array([0.25, 0.25, 0.25, 0.25]) * visit(tree)[1])


def build_emission_matrix(topologies, nleaves, interval_times, theta):
    I = len(interval_times)-1
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

    S = [[None] * I for i in xrange(I)]
    for i in xrange(I):
        for j in xrange(i, I):
            S[i][j] = jukes_cantor(D(i,j))

    res = []
    topology_map = []
    for topo in topologies:
        temp_res = []
        topology_map.append(topo)
        for cols_v in xrange(1 << nleaves*2):
            cols = index_to_cols(cols_v, nleaves)
            row = emission_test(topo, S, cols)
            temp_res.append(row)
        res.append(temp_res)

    # The index of the matrix is [topo, variant (AAA, AAC, etc.)]
    # topology_map gives a mapping from row index to the corresponding tree.
    return matrix(res), topology_map
