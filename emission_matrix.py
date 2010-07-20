from scipy import *
iset = frozenset

def _only(s):
    """At the leaves of a tree, there is always isets of size one.
       This method is used to extract that member."""
    assert len(s) == 1, "All leaves must be of size one"
    for x in s:
        return x

# A column of symbols is converted to an integer by giving each entry
#  2 bits in the result, so: a,b,c,d <=> a + (b << 2) + (c << 4) + (d << 6)
def cols_to_index(cols):
    """Converts a column of up to four symbols so an integer"""
    i = 0
    s = 0
    for v in cols:
        s += v << i
        i += 2
    return s
def index_to_cols(value, n):
    """Converts an integer to a column of n symbols"""
    return map(lambda i: (value >> 2*i) & 3, xrange(n))

# Instead of a matrix we can calculate JC directly
def _jukes_cantor(a, b, dt):
    if a == b:
        return 0.25 + 0.75 * exp(-4*dt)
    else:
        return 0.25 - 0.25 * exp(-4*dt)

def _emission_row(tree, cols, times, theta):
    # m calculates the time in an interval that should be used for further
    # calculation. There is a special case at the end where the time goes to
    # infinity, where we return theta.
    def m(i):
        if i + 1 >= len(times):
            dt = 0.0
            a = 0.0
        else:
            dt = times[i+1] - times[i]
            a = exp(-dt/theta)
        return theta - (dt * a)/(1 - a)

    # Calculates emission probs for a tree.
    # This is done bottom-up, by giving each leaf 1.0 for the symbol it
    # represents.
    # For each subtree in a node, a double sum is needed:
    #  for each symbol, y, in our result:
    #    for each symbol, x, we could be coming from:
    #      add the prob of the child being in x, times the prob of going from
    #      x to y in the time that has passed
    # At the end the node takes the product of the results for each child.
    def visit(t):
        if isinstance(t, iset): # leaf
            res = array([0.0, 0.0, 0.0, 0.0])
            res[cols[_only(t)]] = 1.0
            return 0.0, res
        else: # node
            j, children = t
            m_j = m(j)
            prob_j = ones(4)
            for m_i, prob_i in map(visit, children):
                dt = m_j - m_i
                prob = zeros(4)
                for y in xrange(4):
                    for x in xrange(4):
                        prob[y] += prob_i[x] * _jukes_cantor(x,y,dt)
                prob_j = prob_j * prob
            return m_j, prob_j
    # After calculating the emission probs for the tree, we multiply by
    # a set of weights, currently 1/4 in all entries.
    return sum(array([0.25, 0.25, 0.25, 0.25]) * visit(tree)[1])


def build_emission_matrix(topologies, tmap, nleaves, interval_times, theta):
    I = len(interval_times)
    res = [None] * len(tmap)
    for topo in topologies:
        temp_res = []
        for cols_v in xrange(1 << nleaves*2):
            cols = index_to_cols(cols_v, nleaves)
            row = _emission_row(topo, cols, interval_times, theta)
            temp_res.append(row)
        res[tmap[topo]] = temp_res

    # The index of the matrix is [topo, variant (AAA, AAC, etc.)]
    return matrix(res)
