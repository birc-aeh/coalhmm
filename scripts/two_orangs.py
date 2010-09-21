import sys
import os
sys.path.append("../")
from optimize import logLikelihood, readObservations
from model import build_simple_model, build_epoch_seperated_model
from fasta_parser import readAlignment
from scipy import linspace
from scipy.optimize import fmin

Ne = 25000.0
gen = 20
theta = 2*Ne * gen * 1e-9

C = 1.0 / theta
R = (1.5e-8/gen) / 1.0e-9
tau = 325000e-9

def popsize_to_C(ne):
    return 1.0 / (2.0 * gen * ne * 1e-9)
def C_to_popsize(c):
    return 1.0 / c / (2 * gen * 1e-9)
def recrate_to_R(rec):
    return rec / (gen * 0.1)
def R_to_recrate(r):
    return r * gen * 0.1

def callback(x):
    #print x
    pass

def logL(model, obs, c, r, t):
    if min([c,r,t]) < 0:
        return -1e18
    return logLikelihood(model, obs, [c,c], [r,r], [0,0], [0.0,t])

def optimize(file_in, model, seqnames, init, cb=None):
    obs = readObservations(file_in, seqnames)
    nbps = model.nbreakpoints
    return fmin(lambda x: -logL(model, obs, *x), init, callback=cb)

if __name__ == "__main__":
    folder = "two_orangs"
    _c = int(sys.argv[1])
    variant = 1
    for nstates in [int(sys.argv[2])]:
        noBrPointsPerEpoch = [1, nstates]
        migs = None#[[[1],[0]],None]
        maps = [[0,0]]
        model = build_epoch_seperated_model(2, maps, noBrPointsPerEpoch, migs)
        f = open("%s/nomig%s_%i-%istates.txt" % (folder,_c,variant,nstates), 'w')
        print "Writing to:", f.name
        print >>f, "C\tR\tM\ttau_1"
        for run in xrange(1,26):
            filename = "../../simulated_data/Mig%i_1/Sim%i/seq.unif.fasta" % (_c, run)
            c,r,t = optimize(filename, model, ["'0'","'1'"], (C,R,tau), callback)
            print >>f, '\t'.join(map(repr, [c,r,0.0,t]))
            print "Done with", run
            f.flush()
        f.close()

   # samples = 15
   # Cs = map(popsize_to_C, linspace(0.1*Ne, 2.0*Ne, samples))
   # Rs = map(recrate_to_R, linspace(0.1, 3.0, samples))
   # Ts = linspace(0.1*tau, 2.0*tau ,samples)
   # Ms = linspace(0.1/samples, 0.1, 10*samples)

   # filename = "../../simulated_data/Mig250_1/Sim10/seq.unif.fasta"

   # obs = readObservations(filename, ["'1'", "'2'"])

   # for nstates in [2]:
   #     noBrPointsPerEpoch = [1, nstates, nstates]
   #     migrations = [None,[[1],[0]],None]
   #     maps = [[0,1],[0,0]]
   #     model = build_epoch_seperated_model(2, maps, noBrPointsPerEpoch, migrations)
   #     logFile = open('%s/mle_C_%i.txt' % (folder,nstates),'w')
   #     print >>logFile, '\t'.join(['Recrate', 'Popsize', 'Migration', 'tau', 'logL'])
   #     for t in [tau]:
   #         for c in Cs:
   #             for r in [R]:
   #                 for m in [0.01]:
   #                     print >>logFile, '\t'.join(map(repr,
   #                      [R_to_recrate(r),
   #                          C_to_popsize(c),
   #                          m,
   #                          t,
   #                          logLikelihood(model, obs,
   #                              [c, c, c],
   #                              [r, r, r],
   #                              [m, m, m],
   #                              [0, 0.25*t, t])]))
   #                 logFile.flush()
   #     logFile.close()
