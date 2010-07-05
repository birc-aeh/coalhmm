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
    #print containing_c
    for i in xrange(len(V)):
        comp = V[i]
        for a in comp:
            for b in edges[a]:
                E[i].add(containing_c[b])
                #E[containing_c[b]].add(i)
        E[i].remove(i)

    return V, E
