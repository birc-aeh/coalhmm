from scipy import array, zeros, exp, ones, matrix
from statespace_generator import IM
from coal_time_computer import CoalTimeComputer
iset = frozenset

def _only(s):
    """At the leaves of a tree, there is always isets of size one.
       This method is used to extract that member."""
    assert len(s) == 1, "All leaves must be of size one"
    for x in s:
        return x

# Instead of a matrix we can calculate JC directly
def _jukes_cantor(a, b, dt):
    if a == b:
        return 0.25 + 0.75 * exp(-4.0/3*dt)
    else:
        return 0.25 - 0.25 * exp(-4.0/3*dt)

_to_val = {'A':0,'C':1,'G':2,'T':3}
def _leaf_prob(cols, species):
    symbol = cols[species]
    if symbol in ['N', '-']:
        return ones(4)
    else:
        res = array(zeros(4))
        res[_to_val[symbol]] = 1.0
        return res

def _emission_row(tree, cols, cost):
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
            return 0, _leaf_prob(cols, _only(t))
        else: # node
            j, children = t
            prob_j = ones(4)
            for i, prob_i in map(visit, children):
                prob = zeros(4)
                for y in xrange(4):
                    for x in xrange(4):
                        prob[y] += prob_i[x] * cost[i,j,int(x!=y)]
                prob_j = prob_j * prob
            return j, prob_j
    # After calculating the emission probs for the tree, we multiply by
    # a set of weights, currently 1/4 in all entries.
    return 0.25*sum(visit(tree)[1])


def build_emission_matrix(topologies, tmap, col_map, nleaves,
        interval_times, interval_to_epoch, theta, Qs, G, rates):
    npossible_cols = len(set(col_map.values()))

    def coalTimeSimple(t1, t2):
        dt = t2 - t1
        a = exp(-dt/theta)
        return t1 + theta - (dt * a)/(1 - a)

    nepochs = len(Qs)
    cts = {}
    for e in range(nepochs):
        cts[e] = coalTimeSimple
    for e in range(nepochs):
        if e in cts:
            ct = cts[e]
        else:
            e_start = 0
            for j in range(len(interval_to_epoch)):
                if interval_to_epoch[j] == e:
                    e_start = j
                    break
            e_G = G.G[e]
            rate = rates[interval]
            if e_G.has_migration():
                states, transitions = IM(range(2)).compute_state_space()
                ct = CoalTimeComputer(rate, states, transitions, interval_times[e_start])
                cts[e] = ct

    # m calculates the time in an interval that should be used for further
    # calculation. There is a special case at the end where the time goes to
    # infinity, where we return theta.
    def m(i, times):
        if i + 1 >= len(times):
            return times[-1] + theta
        else:
            t1 = times[i]
            t2 = times[i+1]
            return cts[interval_to_epoch[i]](t1, t2)
    pp_m = [m(i, interval_times) for i in range(len(interval_times))]
    cost = zeros((nepochs,nepochs,2))
    # our cost function can be easily precalculated since we use JC69
    for i in range(nepochs):
        for j in range(i,nepochs):
            for o in [0,1]:
                dt = pp_m[j] - pp_m[i]
                cost[i,j,o] = _jukes_cantor(0,o,dt)
    res = zeros((len(tmap), npossible_cols))
    for topo, topo_i in tmap.iteritems():
        temp_res = []
        for cols, cols_v in col_map.iteritems():
            row = _emission_row(topo, cols, cost)
            res[topo_i, cols_v] += row

    # The index of the matrix is [topo, variant (AAA, AAC, etc.)]
    return matrix(res)
