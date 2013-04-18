# This example is a copy of estimate_example.py, except that it uses zipHMM
# as its HMM implementation. zipHMM needs to preprocess each input sequence
# once, but can then calculate the forward algorithm much faster which is
# what we do when we estimate parameters.

# I assume that zipHMM is either installed or that there is a local copy of
# pyZipHMM.py and libpyZipHMM.so in the same folder.


from coalhmm.optimize import logL_multiseq
# logL_multiseq is made to work with other HMM libraries, we just need to
# provide two functions. One that maps from scipy matrices to whatever the new
# HMM library uses and one that evaluates the forward algorithm given the
# transformed matrices and a single sequence.
from pyZipHMM import *

# zipHMM has its own Matrix class, so we copy from the scipy matrix in to one
# of those.
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
# in zipHMM the observations know how to do the forward algorithm.
def zipHMM_single_logL(pi, T, E, obs):
    return obs.forward(pi, T, E)

# logLikelihood is changed to use our two new functions.
def logLikelihood(model, all_obs, c,r,m,t):
    return logL_multiseq(model, all_obs, FIXED_COL_MAP, c,r,m,t,
            prepare_matrices=zipHMM_prepare_matrices,
            single_logL=zipHMM_single_logL)


current_model = 'I'

intervals_per_epoch = 10

from coalhmm.model import build_epoch_separated_model
if current_model == 'I':
    modelI = build_epoch_separated_model(
            2,
            [[0,0]],
            [1,intervals_per_epoch])
elif current_model == 'IM':
    M_for_both = [[1],[0]]
    modelIM = build_epoch_separated_model(
            2,
            [[0,1], [0,0]],
            [1,intervals_per_epoch,intervals_per_epoch],
            [None, M_for_both, None])

names = ["hg18", "pantro2"]

from coalhmm.fasta_parser import readAlignment

from scipy import zeros, int32
def read_observations(filename, names):
    '''Reads in a single fasta file and converts it to a series of
    observations.  The observations assume JC69, so we only care if there is a
    difference or not (or a N).'''
    alignments = readAlignment(filename)
    srcs = [alignments[name] for name in names]
    # alignments = SeqIO.to_dict(SeqIO.parse(open(filename),'fasta'))
    # srcs = [alignments[name].seq for name in names]

    clean = set(['A', 'C', 'G', 'T'])
    A = srcs[0]
    B = srcs[1]
    assert len(A) == len(B)
    L = len(A)
    # int16/int8 can save some memory when using larger datasets
    obs = zeros(L, dtype=int32) 
    for i in xrange(L):
        s1,s2 = A[i], B[i]
        if s1 not in clean or s2 not in clean:
            v = 2
        elif s1 == s2:
            v = 0
        else:
            v = 1
        obs[i] = v
    return obs

# Unfortunately zipHMM needs to read sequences from a file with a series of
# observations(numbers) separated by spaces.
# So after reading in our observations we write them back out again.
import tempfile
import os
all_obs = []
total_L = 0
for seq in ["example_data.fa"]:
    obs = read_observations(seq, names)
    total_L += len(obs)
    fd, foutname = tempfile.mkstemp()
    fout = os.fdopen(fd, 'w', 64*1024)
    for o in obs:
        print >>fout, o,
    fout.close()
    if current_model == 'I':
        # when constructing the sequence we need to tell it the number of
        # possible observations and the number of states in the HMM (the size
        # of the transition matrix).
        f = Forwarder(seqFilename = foutname, nStates = len(modelI.tree_map), nObservables = 3)
    elif current_model == 'IM':
        f = Forwarder(seqFilename = foutname, nStates = len(modelIM.tree_map), nObservables = 3)
    all_obs.append(f)
    os.remove(foutname)
print "%.2fMbp of data in total" % (total_L/1e6)

g = 20   # generation time
u = 1e-9 # mutation rate
# initial values for recombination rate, population size/coalescence rate and
# migration rate:
i_r = 0.8
i_N = 20e3
i_c = 1.0/(2*g*u*i_N)
i_m = 0.1

init_I =   (5.0e6*u, i_c, i_r)
init_IM =  (4.0e6*u, 5.0e6*u, i_m/2*i_c, i_c, i_r)

FIXED_COL_MAP = { ('A','A'): 0, ('A','C'): 1, ('N','N'): 2 }

from scipy.optimize import fmin
def optimize_f(f_logL, init):
    return fmin(lambda x: -f_logL(*x), init, full_output=True)

def estimate_I(model, all_obs, T, C, R):
    def logL_all(t, c, r):
        if min([t,c,r]) <= 0:
            return -1e18
        return logLikelihood(model, all_obs, [c]*2, [r]*2, [0]*2, [0.0,t])
    est, L, _, _, _ = optimize_f(logL_all, (T,C,R))
    return (L, list(est))

def estimate_IM(model, all_obs, T1, T2, M, C, R):
    def logL_all(t1, t2, m, c, r):
        if min([t1,t2,m,c,r]) <= 0 or t2 <= t1:
            return -1e18
        return logLikelihood(model, all_obs, [c]*3, [r]*3, [0.0,m,0.0], [0.0,t1,t2])
    est, L, _, _, _ = optimize_f(logL_all, (T1,T2,M,C,R))
    return (L, list(est))

if current_model == 'I':
    L, est = estimate_I(modelI, all_obs, *init_I)
elif current_model == 'IM':
    L, est = estimate_IM(modelIM, all_obs, *init_IM)

print 'Final results:'
print "\t".join(map(str, [L] + est))
