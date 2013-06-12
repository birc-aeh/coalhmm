import os, sys
from optparse import OptionParser

from pyZipHMM import Matrix, Forwarder

from time import time
__t0 = time()
def now():
    return time() - __t0

import coalhmm
from coalhmm.model import build_epoch_seperated_model
from coalhmm.optimize import logL_multiseq
from scipy.optimize import fmin

def optimize_f(model, obs, f_logL, init):
    return fmin(lambda x: -f_logL(model, obs, *x), init, full_output=True)

FIXED_COL_MAP = {
            ('A','A'): 0,
            ('A','C'): 1,
            ('N','N'): 2
            }

def zipHMM_prepare_matrices(mpi, mT, mE):
    pi = Matrix(mpi.shape[0], 1)
    for i in xrange(mpi.shape[0]):
        pi[0, i] = mpi[i]
    T  = Matrix(mT.shape[0], mT.shape[1])
    for i in xrange(mT.shape[0]):
        for j in xrange(mT.shape[1]):
            T[i, j] = mT[i, j]
    E  = Matrix(mE.shape[0], mE.shape[1])
    for i in xrange(mE.shape[0]):
        for j in xrange(mE.shape[1]):
            E[i, j] = mE[i, j]
    return pi, T, E

def zipHMM_single_logL(pi, T, E, obs):
    return obs.forward(pi, T, E)

def logLikelihood(model, all_obs, c,r,m,t):
    return logL_multiseq(model, all_obs, FIXED_COL_MAP, c,r,m,t,
            prepare_matrices=zipHMM_prepare_matrices,
            single_logL=zipHMM_single_logL)


def estimate_I(model, obs, T, C, R, outfile):
    def logL_all(model, all_obs, t, c, r):
        if min([t,c,r]) <= 0:
            return -1e18
        res = logLikelihood(model, all_obs, [c]*2, [r]*2, [0.0]*2, [0.0,t])
        print >>outfile, '\t'.join(map(str, [now(), t, c, r, res]))
        return res
    print >>outfile, '\t'.join(map(str, ["time", "t", "c", "r", "logL"]))
    est, L, _, _, _ = optimize_f(model, obs, logL_all, (T,C,R))
    #       logL   t,c,r
    return (L,     est)

def log_unfinished_line(s):
    print s,
    sys.stdout.flush()
def log_finished_line(s):
    print s

def main():
    usage="""%prog [options] <forwarder dirs>

This program estimates the parameters of an isolation model with two species
and uniform coalescence/recombination rate."""


    parser = OptionParser(usage=usage, version="%prog 1.0")

    parser.add_option("-o", "--out",
                      dest="outfile",
                      type="string",
                      default="/dev/stdout",
                      help="Output file for the estimate (/dev/stdout)")
    parser.add_option("--tmpfile",
                      dest="tmpfile",
                      type="string",
                      default="/dev/null",
                      help="Log for all points estimated in the optimization (/dev/null)")
    optimized_params = [
            ('coaltime', 'coalescence time', 1e6),
            ('Ne', 'effective population size', 20e3),
            ('recomb', 'recombination rate', 0.1),
            ]
    for (cname, desc, default) in optimized_params:
        parser.add_option("--%s" % cname,
                          dest=cname,
                          type="float",
                          default=default,
                          help="Initial guess at the %s (%g)" % (desc, default))
    fixed_params = [
            ('mu', 'mutation rate', 1e-9),
            ('g', 'generation time', 20),
            ]
    for (cname, desc, default) in fixed_params:
        parser.add_option("--%s" % cname,
                          dest=cname,
                          type="float",
                          default=default,
                          help="Value of the %s (%g)" % (desc, default))
    parser.add_option("--intervals",
                      dest="intervals",
                      type="int",
                      default=10,
                      help="Number of sub intervals used to discretize the time (10)")
    parser.add_option("--header",
                      dest="include_header",
                      action="store_true",
                      default=False,
                      help="Include a header on the output")
    parser.add_option("-v", "--verbose",
                      dest="verbose",
                      action="store_true",
                      default=False,
                      help="Print some stuff")

    (options, args) = parser.parse_args()
    if len(args) < 1:
        parser.error("Needs at least one preprocessed sequence to work on")

    if not options.verbose:
        log = lambda s: None
        logu = lambda s: None
    else:
        logu = log_unfinished_line
        log = log_finished_line

    logu("Loading forwarders...")
    forwarders = [Forwarder(dir) for dir in args]
    log("done")

    logu("Constructing model...")
    intervals = options.intervals
    modelI = build_epoch_seperated_model(2, [[0,0]], [1,intervals])
    log("done")


    mu = options.mu
    g = options.g
    T = options.coaltime * mu
    C = 1.0/(g*mu*2*options.Ne)
    R = options.recomb
    with open(options.tmpfile, 'w') as tmpfile:
        L, est = estimate_I(modelI, forwarders, T, C, R, outfile=tmpfile)
    vals = "\t".join(map(str,est))
    with open(options.outfile, 'w') as outfile:
        if options.include_header:
            print >>outfile, 'logL\tTcoal\tC\tR'
        print >>outfile, "%f\t%s" % (L,vals)

if __name__ == "__main__":
    main()
