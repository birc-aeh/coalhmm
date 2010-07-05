from sets import ImmutableSet as iset

class Tarjan:
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


def build_scc(states, edges):
    t = Tarjan(len(states), edges)
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

    return V, E

class SCCGraph:
    def __init__(self, states, edges):
        self.states_rev = {}
        for k, v in states.iteritems():
            self.states_rev[v] = k

        edges2 = [[] for i in xrange(len(states))]
        for (a,t,b) in edges:
            edges2[a].append(b)
        self.V, self.E = build_scc(states, edges2)

    def add_transitive_edges(self):
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

    def initial(self):
        return self.states_rev[0]
    def state(self, v):
        return self.states_rev[self.V[v][0]]
    def project_state(self, s, side):
        return iset([x[side] for x in s if x[side] != iset()])
    def projected(self, v, side):
        return self.project_state(self.state(v), side)

