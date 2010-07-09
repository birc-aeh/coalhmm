from sets import ImmutableSet as iset

# A standard implementation of Tarjans algorithm.
# The class is just to keep track of some state instead of having a lot of
# params.
class _Tarjan:
    def __init__(self, N, edges):
        self.S = []
        self.SCC = []
        self.index = [-1] * N
        self.low = [-1] * N
        self.max_index = 0
        self.edges = edges
    def tarjan(self, a):
        index = self.index
        low = self.low
        S = self.S
        if low[a] != -1:
            return
        index[a] = low[a] = self.max_index
        self.max_index = self.max_index + 1
        S.append(a)
        for b in self.edges[a]:
            if index[b] == -1:
                self.tarjan(b)
                low[a] = min(low[a],low[b])
            elif S.count(b) == 1:
                low[a] = min(low[a],index[b])
        if low[a] == index[a]:
            n = S.pop()
            component = [n]
            while n != a and len(S) > 0:
                n = S.pop()
                component.append(n)
            self.SCC.append(component)


def _build_scc(states, edges):
    t = _Tarjan(len(states), edges)
    for v in xrange(len(states)):
        t.tarjan(v)
    V = t.SCC
    E = [set() for c in V]
    containing_c = [-1] * len(states)
    for i in xrange(len(V)):
        for v in V[i]:
            containing_c[v] = i
    for i in xrange(len(V)):
        comp = V[i]
        for a in comp:
            for b in edges[a]:
                E[i].add(containing_c[b])
        E[i].remove(i)

    return V, E, containing_c

class SCCGraph:
    def merge(Ga, Gb):
        # All states in Ga and Gb are added, with a new number.
        # Ga receives the same numbers, while Gb states add len(Ga's states)
        offset = len(Ga.original_states)
        new_states = dict(Ga.original_states)
        new_states_rev = dict(Ga.states_rev)
        for k, v in Gb.original_states.iteritems():
            new_states[k] = v + offset
            new_states_rev[v + offset] = k
        new_edges = Ga.original_edges[:]
        for (a,t,b) in Gb.original_edges:
            new_edges.append((a+offset, t, b+offset))
        new_V = Ga.V[:]
        for c in Gb.V:
            new_V.append([x+offset for x in c])
        offset = len(Ga.V) # now we are using component offsets
        new_E = Ga.E[:]
        for edges in Gb.E:
            new_E.append(set([e+offset for e in edges]))
        new_cc = Ga.containing_c + [x + offset for x in Gb.containing_c]

        new_G = SCCGraph(None, None)
        new_G.original_states = new_states
        new_G.original_edges = new_edges
        new_G.states_rev = new_states_rev
        new_G.V = [sorted(x) for x in new_V]
        new_G.E = new_E
        new_G.containing_c = new_cc
        new_G.epochs = Ga.epochs + Gb.epochs
        return new_G

    def __init__(self, states, edges, epoch=0):
        if states == None and edges == None:
            return
        self.original_edges = edges
        self.original_states = states
        self.states_rev = {}
        for k, v in states.iteritems():
            self.states_rev[v] = k

        edges2 = [[] for i in xrange(len(states))]
        for (a,t,b) in edges:
            edges2[a].append(b)
        self.V, self.E, self.containing_c = _build_scc(states, edges2)
        self.epochs = [epoch] * len(self.V)

    def originalGraph(self):
        return self.original_states, self.original_edges

    def add_transitive_edges(self):
        """Give each component extra edges to everything that is reachable"""
        def transitive_edges(a, S, newE):
            for v in S:
                newE[v].add(a)
            S.append(a)
            for b in self.E[a]:
                transitive_edges(b, S, newE)
            S.pop()
        E2 = [set() for x in self.E]
        transitive_edges(len(self.V)-1, [], E2)
        self.E = E2
    def add_component_edge(self, c1, c2):
        assert isinstance(c1, iset) and isinstance(c2, iset) and c1 != c2
        res = self.original_states[c1], self.original_states[c2]
        c1, c2 = self.find_component(c1), self.find_component(c2)
        self.E[c1].add(c2)
        return res

    def find_component(self, state):
        return self.containing_c[self.original_states[state]]

    def initial(self):
        """Returns the initial state, where all species are separate"""
        return self.states_rev[0]
    def all_states(self, v):
        """Returns all the states in a component"""
        return [x for x in self.V[v]]
    def state(self, v):
        """Returns a representative state for a component"""
        return self.states_rev[self.V[v][0]]
    def project_state(self, s, side):
        """Extracts either the left or right side from a state"""
        return iset([x[1][side] for x in s if x[1][side] != iset()])
    def projected(self, v, side):
        """Extracts either the left or right side from a representative of the
        component"""
        return self.project_state(self.state(v), side)

