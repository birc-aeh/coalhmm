import os, sys
from optparse import OptionParser

from pyZipHMM import Matrix, posteriorDecoding
from scipy import identity

import coalhmm
from coalhmm.model import build_epoch_seperated_model
from coalhmm.optimize import default_bps

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

def get_model_matrices(model, C, R, T):

    c = [C] * 2
    r = [R] * 2
    t = [0.0, T]
    
    noBrPointsPerEpoch = model.nbreakpoints
    nleaves = model.nleaves
    all_time_breakpoints, time_breakpoints = default_bps(model, c, r, t)

    M = []
    for e in xrange(len(noBrPointsPerEpoch)):
        newM = identity(nleaves)
        newM[:] = 0
        M.append(newM)

    pi, T, E = model.run(r, c, time_breakpoints, M, col_map=FIXED_COL_MAP)
    return pi, T, E

def log_unfinished_line(s):
    print s,
    sys.stdout.flush()
def log_finished_line(s):
    print s

def main():
    usage="""%prog [options] <forwarder dirs>

This program calculates the posterior state probabilities of an isolation
model with two species and uniform coalescence/recombination rate."""


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
            ('T', 'split time'),
            ('C', 'Coalescence rate'),
            ('R', 'recombination rate'),
            ]
    for (cname, desc) in optimized_params:
        parser.add_option("-%s" % cname,
                          dest=cname,
                          type="float",
                          help="Model parameter %s" % desc)
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
                      help="Print help")

    (options, args) = parser.parse_args()
    if len(args) < 1:
        parser.error("Needs at least one preprocessed sequence to work on")

    if not options.verbose:
        log = lambda s: None
        logu = lambda s: None
    else:
        logu = log_unfinished_line
        log = log_finished_line

    logu("Constructing model...")
    intervals = options.intervals
    model = build_epoch_seperated_model(2, [[0,0]], [1,intervals])
    log("done")


    T = options.T 
    C = options.C
    R = options.R
    
    if T is None:
        print 'You must specify the split time.'
        sys.exit(2)
    if C is None:
        print 'You must specify the coalescence rate.'
        sys.exit(2)
    if R is None:
        print 'You must specify the recombination rate.'
        sys.exit(2)
    

    logu("Preparing the HMM...")    
    mpi, mT, mE = get_model_matrices(model, C, R, T)
    pi,   T,  E = zipHMM_prepare_matrices(mpi, mT, mE)
    log("done")


    for ziphmmdir in args:    
        logu("Computing posterior table for %s..."%ziphmmdir)
        pdPath, pdTable = posteriorDecoding(ziphmmdir, pi, T, E)
        log("done")
    
    
    

if __name__ == "__main__":
    main()
