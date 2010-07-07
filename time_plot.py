from scipy.linalg import expm
from scipy import *
from pylab import *
from scc import SCCGraph

def prettify_state(s):
    """Convert a coal system state to something nicer.
    
    example:
      iset([(iset([3]), iset([3])),
            (iset([1]), iset([1])),
            (iset([2]), iset([2]))])
    to:
      {3, 1, 2}, {3, 1, 2}
    """
    def f(s, side, d):
        if d == 0:
            tmp = [f(sub, side, d+1) for sub in s if len(s) > 0]
            return ", ".join([x for x in tmp if x.strip() != ""])
        elif d == 1:
            return "".join(sorted([str(x) for x in s[side] if len(s[side]) > 0]))
    return "{" + f(s, 0, 0) + "}, {" + f(s, 1, 0) + "}"

def time_plot(G, Q, T):
    times = linspace(0, T)
    data = array([expm(Q*t)[0,:] for t in times])
    for v in xrange(len(G.V)):
        plot(times, sum(data[:,G.all_states(v)], axis=1), label=prettify_state(G.state(v)))
        legend()
        savefig("timeplots/comp_"+str(v)+".png")
        clf()

