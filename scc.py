iset = frozenset

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
        if i in E[i]:
            E[i].remove(i)

    return V, E, containing_c

class EpochSeperatedSCCGraph:
    def __init__(self):
        self.G = []
        self.E = {}

    def addSubGraph(self, sg):
        self.G.append(sg)

    def getEpochSizes(self):
        return [len(x.original_states) for x in self.G]

    def originalGraph(self, epoch):
        return self.G[epoch].originalGraph()

    def add_component_edge(self, e1, s1, e2, s2):
        assert e1 != e2
        res = self.G[e1].original_states[s1], self.G[e2].original_states[s2]
        c1, c2 = self.G[e1].find_component(s1), self.G[e2].find_component(s2)
        self.E[(e1,c1)] = (e2,c2)

        return res

    def initial(self, epoch):
        """Returns the initial state, where all species are separate"""
        return self.G[epoch].initial()
    def all_states(self, epoch, v):
        """Returns all the states in a component"""
        return self.G[epoch].all_states(v)
    def state(self, epoch, v):
        """Returns a representative state for a component"""
        return self.G[epoch].state(v)
    def project_state(self, epoch, s, side):
        """Extracts either the left or right side from a state"""
        return self.G[epoch].project_state(s, side)
    def projected(self, epoch, v, side):
        """Extracts either the left or right side from a representative of the
        component"""
        return self.G[epoch].projected(v, side)

    def all_paths(self):
        allE = dict([(k, set([v])) for k, v in self.E.iteritems()])
        for e in xrange(len(self.G)):
            Ge = self.G[e]
            for c, links in enumerate(Ge.E):
                for link in links:
                    allE.setdefault((e,c), set()).add((e,link))
        S = []
        def dfs(a):
            S.append(a)
            if a not in allE:
                yield S
            else:
                for b in allE[a]:
                    for res in dfs(b):
                        yield res
            S.pop()
        return dfs((0,0))

class SCCGraph:
    def edges_at_component_level(self):
        """Write a dot format graph of this graph to f."""
        for c, links in enumerate(self.E):
            for link in links:
                yield c, link

    def __init__(self, states, edges, epoch=0):
        if states == None and edges == None:
            return
        self.original_edges = edges
        self.original_states = states
        self.states_rev = {}
        for k, v in states.iteritems():
            self.states_rev[v] = k

        edges2 = [[] for i in xrange(len(states))]
        for (a,t,pop_a,pop_b,b) in edges:
            edges2[a].append(b)
        self.V, self.E, self.containing_c = _build_scc(states, edges2)
        self.projection_cache = {}
        self.state_tuples = [tuple(self.V[v]) for v in xrange(len(self.V))]

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

    def find_component(self, state):
        return self.containing_c[self.original_states[state]]

    def initial(self):
        """Returns the initial state, where all species are separate"""
        return self.states_rev[0]
    def all_states(self, v):
        """Returns all the states in a component"""
        return self.state_tuples[v]
    def state(self, v):
        """Returns a representative state for a component"""
        return self.states_rev[self.V[v][0]]
    def project_state(self, s, side):
        """Extracts either the left or right side from a state"""
        if (s, side) not in self.projection_cache:
            res = iset([x[1][side] for x in s if x[1][side] != iset()])
            self.projection_cache[(s,side)] = res
            return res
        return self.projection_cache[(s,side)]
    def projected(self, v, side):
        """Extracts either the left or right side from a representative of the
        component"""
        return self.project_state(self.state(v), side)

