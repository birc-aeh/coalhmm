import sys
import os
sys.path.append("../")
from optimize import logLikelihood, readObservations
from scipy.optimize import fmin
from model import build_simple_model, build_epoch_seperated_model
from fasta_parser import readAlignment
from itertools import izip

Ne = 30000.0
gen = 20
theta = 2*Ne * gen * 1e-9

C = 1.0 / theta
R = (1.5e-8/gen) / 1.0e-9
tau = 500000e-9

def popsize_to_C(ne):
    return 1.0 / (2.0 * gen * ne * 1e-9)
def C_to_popsize(c):
    return 1.0 / c / (2 * gen * 1e-9)
def recrate_to_R(rec):
    return rec / (gen * 0.1)
def R_to_recrate(r):
    return r * gen * 0.1

def logL(model, obs, c, r, t):
    if c < 0 or r < 0 or t < 0:
        return -1e18
    return logLikelihood(model, obs, [c, c], [r, r], [0, 0], [0.0, t])

def optimize(file_in, model, seqnames, init):
    obs = readObservations(file_in, seqnames)
    nbps = model.nbreakpoints
    return fmin(lambda x: -logL(model, obs, *x), init)

if __name__ == "__main__":
    print "From params:", 1e-9*(0.5e6 + 2*Ne*gen)
    model = build_epoch_seperated_model(2, [[0,0]], [1,4])
    f = open("res2.txt", "w")
    print >>f, "counted\tC\tR\ttau_1"
    for i in ["01", "02", "03", "04", "05", "06",
              "07", "08", "09", "10", "11"]:
        filename = "20101005_sim%s" % i
    #for i in [1]:
    #    filename = "../simulated_20101005/sim_4_0.txt"
        alignments = readAlignment(filename)
        seqnames = ["'1'", "'2'"]
        ca = alignments[seqnames[0]]
        cb = alignments[seqnames[1]]
        diffs = 0
        for a, b in izip(ca, cb):
            if a != b:
                diffs += 1
        print "Counted:", 1.0 * diffs / len(ca)
        c,r,t = optimize(filename, model, seqnames, (C,R,tau))
        print "From model estimate:", (t + 1/c), c, r, t
        print >>f, "\t".join(map(str, [1.0*diffs / len(ca), c, r, t]))
        f.flush()
    f.close()
