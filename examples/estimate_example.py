import sys
import os

import coalhmm
from coalhmm.model import build_epoch_seperated_model
from coalhmm.optimize import logL_multiseq
from coalhmm.fasta_parser import readAlignment
from scipy import *
from scipy.optimize import fmin


# If zipHMM is installed the whole thing will likely run much faster.
# from pyZipHMM import *
# def zipHMM_prepare_matrices(mpi, mT, mE):
#     pi = Matrix(mpi.shape[0], 1)
#     for i in xrange(mpi.shape[0]):
#         pi[0, i] = mpi[i]
#     T  = Matrix(mT.shape[0], mT.shape[1])
#     for i in xrange(mT.shape[0]):
#         for j in xrange(mT.shape[1]):
#             T[i, j] = mT[i, j]
#     E  = Matrix(mE.shape[0], mE.shape[1])
#     for i in xrange(mE.shape[0]):
#         for j in xrange(mE.shape[1]):
#             E[i, j] = mE[i, j]
#     return pi, T, E
# 
# def zipHMM_single_logL(pi, T, E, obs):
#     return obs.forward(pi, T, E)

# When calculating the model needs a 'column map'. It maps from a column of
# symbols to a value in the sequence of observations.
# In this example we assume JC69 so we have three values, one for two equal
# symbols, one for two different symbols and one for unknown data.
# An inverse of the map is made so it must be a 1-to-1 mapping.
FIXED_COL_MAP = {
            ('A','A'): 0,
            ('A','C'): 1,
            ('N','N'): 2
            }

def logLikelihood(model, all_obs, c,r,m,t):
    return logL_multiseq(model, all_obs, FIXED_COL_MAP, c,r,m,t)
            #prepare_matrices=zipHMM_prepare_matrices,
            #single_logL=zipHMM_single_logL)

def optimize_f(model, obs, f_logL, init):
    return fmin(lambda x: -f_logL(model, obs, *x), init, full_output=True)

# Here are the functions doing the actual estimating of parameters.
# It works by giving scipy.optimize.fmin a function that evaluates the logL
# given the current parameters. The logging is implemented in an odd way
# that is probably not necessary for your environment.
def estimate_I(model, obs, T, C, R, outfile="/dev/null"):
    def logL_all(model, all_obs, t, c, r):
        if min([t,c,r]) <= 0:
            return -1e18
        res = logLikelihood(model, all_obs, [c]*2, [r]*2, [0]*2, [0.0,t])
        os.system("echo '%s' >> %s" % ('\t'.join(map(str, [t, c, r, res])), outfile))
        return res
    os.system("echo '%s' > %s" % ('\t'.join(map(str, ["t", "c", "r", "logL"])), outfile))
    est, L, _, _, _ = optimize_f(model, obs, logL_all, (T,C,R))
    #       logL.   estimates
    return (L,      list(est))

def estimate_IM(model, obs, T1, T2, M, C, R, outfile="/dev/null"):
    def logL_all(model, all_obs, t1, t2, m, c, r):
        if min([t1,t2,m,c,r]) <= 0 or t2 <= t1:
            return -1e18
        res = logLikelihood(model, all_obs, [c]*3, [r]*3, [0.0,m,0.0], [0.0,t1,t2])
        os.system("echo '%s' >> %s" % ('\t'.join(map(str, [t1, t2, m, c, r, res])), outfile))
        return res
    os.system("echo '%s' > %s" % ('\t'.join(map(str, ["t1", "t2", 'm', "c", "r", "logL"])), outfile))
    est, L, _, _, _ = optimize_f(model, obs, logL_all, (T1,T2,M,C,R))
    #       logL.   estimates
    return (L,      list(est))

def estimate_ILS(model, obs, T1, T2, C, R, outfile="/dev/null"):
    def logL_all(model, all_obs, t1, t2, c, r):
        if min([t1,t2,c,r]) <= 0 or t2 <= t1:
            return -1e18
        res = 0
        for obs, colmap in all_obs:
            res += logLikelihood(model, obs, colmap, [c]*3, [r]*3, [0.0]*3, [0.0,t1,t2])
        os.system("echo '%s' >> %s" % ('\t'.join(map(str, [t1,t2,c, r, res])), outfile))
        return res
    os.system("echo '%s' > %s" % ('\t'.join(map(str, ["t1","t2", "c", "r", "logL"])), outfile))
    est, L, _, _, _ = optimize_f(model, obs, logL_all, (T1,T2,C,R))
    #       logL.   estimates
    return (L, list(est))



def read_observations(filename, names):
    '''Reads in a single fasta file and converts it to a series of
    observations.  The observations assume JC69, so we only care if there is a
    difference or not (or a N).'''
    alignments = readAlignment(filename)
    srcs = [alignments[name] for name in names]

    clean = set(['A', 'C', 'G', 'T'])
    A = srcs[0]
    B = srcs[1]
    assert len(A) == len(B)
    L = len(A)
    # int16 can save some memory when using larger datasets
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


# generation time and mutation rate
g = 20
u = 1e-9

# initial values for recombination rate, population size/coalescence rate and
# migration rate.
i_r = 0.8
i_N = 20e3
i_c = 1.0/(2*g*u*i_N)
i_m = 0.1

estates = 10

current_model = 'I'

name = "%s.test" % (current_model)

# This is just one small chunks of an alignment and won't give very meaningful
# results. Try getting 10Mbp or so per run for simple models.
seqs = ["example_data.fa"]

if current_model == 'I':
    # We can build a simple isolation model with two epochs, in the first epoch
    # everything is separate and at the first epoch transition we map both
    # species to one.
    modelI = build_epoch_seperated_model(
            2,
            [[0,0]],
            [1,estates])
    nstates = len(modelI.tree_map)
    names = ["bonobo", "pantro2"]
    init_I =   (0.5e6*u, i_c, i_r)
elif current_model == 'IM':
    # Creating an isolation-with-migration model is similar to the plain
    # isolation model. We add an extra epoch (maintaining the two separate
    # species) and allow migration from 0 to 1 and from 1 to 0 in it.
    M_for_both = [[1],[0]]
    modelIM = build_epoch_seperated_model(
            2,
            [[0,1], [0,0]],
            [1,estates,estates],
            [None, M_for_both, None])
    nstates = len(modelIM.tree_map)
    names = ["hg18", "pantro2"]
    init_IM =  (4.0e6*u, 5.0e6*u, i_m/2*i_c, i_c, i_r)
elif current_model == 'ILS':
    # Creating an isolation model with three species is done by merging species
    # pairwise.
    # Note: will not work without a modified read_observations method. The one
    # in coalhmm.optimize can be used.
    modelILS = build_epoch_seperated_model(
            3,
            [[0,0,1], [0,0]],
            [1,estates,estates])
    nstates = len(modelILS.tree_map)
    names = ["bonobo", "pantro2", "hg18"]
    init_ILS = (1.0e6*u, 4.5e6*u, i_c, i_r)

all_obs = []
total_L = 0
for seq in seqs:
    obs = read_observations(seq, names)
    total_L += len(obs)
    all_obs.append(obs)

print "%.2fMbp of data in total" % (total_L/1e6)

if current_model == 'I':
    L, est = estimate_I(modelI,     all_obs, *init_I,   outfile="/dev/stdout")
elif current_model == 'IM':
    L, est = estimate_IM(modelIM,   all_obs, *init_IM,  outfile="/dev/stdout")
elif current_model == 'ILS':
    L, est = estimate_ILS(modelILS, all_obs, *init_ILS, outfile="/dev/stdout")
else:
    assert False


print 'Final results:'
print "\t".join(map(str, [L] + est))
