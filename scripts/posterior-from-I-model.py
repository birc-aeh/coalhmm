import os, os.path, sys
import errno
from optparse import OptionParser

from pyZipHMM import Matrix, posteriorDecoding
from scipy import identity
from numpy import array

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
    return pi, T, E, time_breakpoints[1]

def print_matrix(matrix, seqname, first_pos, fout):
    no_states = matrix.getHeight()
    seq_len = matrix.getWidth()

    pos = first_pos
    for line in xrange(seq_len):
        fout.write('%s\t%d' % (seqname, pos))
        pos += 1        
        for state in xrange(no_states):
            fout.write('\t%f' % matrix[state,line])
        fout.write('\n')

def main():
    usage="""%prog [options] <forwarder dir>

This program calculates the posterior state probabilities of an isolation
model with two species and uniform coalescence/recombination rate."""


    parser = OptionParser(usage=usage, version="%prog 1.0")

    parser.add_option("-o", "--out",
                      dest="outfile",
                      type="string",
                      default="/dev/stdout",
                      help="Output file for the estimate (/dev/stdout)")

    parser.add_option("-n", "--seq-name",
                      dest="seq_name",
                      type="string",
                      default="sequence",
                      help="Name of the sequence. Used in the output for tabix indexing")
    parser.add_option("-p", "--first-pos",
                      dest="first_pos",
                      type="int",
                      default=0,
                      help="Position in the sequence where the first element is."
                        "Used in the output for tabix indexing. Default 1.")

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
    parser.add_option("--no-header",
                      dest="no_header",
                      action="store_true",
                      default=False,
                      help="Do not include the header in the output")

    (options, args) = parser.parse_args()
    if len(args) != 1:
        parser.error("Needs a preprocessed sequence to work on")

    split_time = options.T 
    coal_rate  = options.C
    rec_rate   = options.R
    
    if split_time is None:
        print 'You must specify the split time.'
        sys.exit(2)
    if coal_rate is None:
        print 'You must specify the coalescence rate.'
        sys.exit(2)
    if rec_rate is None:
        print 'You must specify the recombination rate.'
        sys.exit(2)
    
    intervals = options.intervals
    model = build_epoch_seperated_model(2, [[0,0]], [1,intervals])

    mpi, mT, mE, time_breaks = get_model_matrices(model, coal_rate, rec_rate, split_time)
    pi,   T,  E = zipHMM_prepare_matrices(mpi, mT, mE)

    with open(options.outfile,'w') as fout:
    
        if not options.no_header:
            print >> fout, '## Posterior probabilities for isolation model.'
            print >> fout, '# intervals =', intervals
            print >> fout, '# T =', split_time
            print >> fout, '# C =', coal_rate
            print >> fout, '# R =', rec_rate
            print >> fout, '# time_breaks =', ' '.join(map(str, time_breaks))
    
        try:
            for ziphmmdir in args:    
                seqfile = os.path.join(ziphmmdir, 'original_sequence')
                _, pdTable = posteriorDecoding(seqfile, pi, T, E)
                print_matrix(pdTable, options.seq_name, options.first_pos, fout)
                
        except IOError as e:
            if e.errno == errno.EPIPE:
                sys.exit(0) # the pipe died, this is probably because it is in
                            # a shell pipe where we have stopped reading
            else:
                raise e

    
    
    

if __name__ == "__main__":
    main()
