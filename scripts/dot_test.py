import sys
import os
sys.path.append("../")
from model import build_simple_model, build_epoch_seperated_scc
from scipy import linspace

if __name__ == "__main__":
    colors = ["red", "blue", "green"]
    prog, g = build_epoch_seperated_scc(2, [[0,0]], [[[1],[0]],None])
    only_cross_epoch_edges = False
    n = lambda e, x: "e"+str(e)+"_c"+str(x)
    #L = lambda e, c, parts: ", ".join(map(str, parts))
    #L = lambda e, c, parts: "n="+str(len(parts))
    L = lambda e, c, parts: n(e, c)
    print "digraph {"
    for e in xrange(2):
        print "subgraph g_"+str(e),"{"
        for c, parts in enumerate(g.G[e].V):
            print n(e, c), '[color="'+colors[e]+'",label="'+L(e,c,parts)+'"];'
        if not only_cross_epoch_edges:
            for a, b in g.G[e].edges_at_component_level():
                print n(e, a), "->", n(e, b)
        print "}"
    for (e1,c1), (e2,c2) in g.E.iteritems():
        print n(e1, c1), "->", n(e2,c2), '[color="red"]'
    print "}"
