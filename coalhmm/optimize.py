from scipy import *
from scipy.stats import expon
from model import build_simple_model, build_epoch_seperated_model
from fasta_parser import readAlignment
from mini_hmm import *
#from pyhmmlib import *
from scipy.optimize import fmin
import sys
import os.path
import time
from itertools import izip

def readObservations(filename, seq_names):
    alignments = readAlignment(filename)
    srcs = [(i, alignments[name]) for i, name in enumerate(seq_names)]

    legal = set(['A', 'C', 'G', 'T', 'N', '-'])
    col_map = dict()
    first = alignments[seq_names[0]]
    #obs = Sequence(len(first))
    obs = array([0]*len(first), dtype=int16)
    tmp = [0 for n in seq_names]
    for i in xrange(len(first)):
        for j, src in srcs:
            s = src[i]
            tmp[j] = s in legal and s or 'N'
        col = tuple(tmp)
        v = col_map.setdefault(col, len(col_map))
        obs[i] = v
    return obs, col_map

def copyTable(dst, src):
   for i in xrange(src.shape[0]):
       for j in xrange(src.shape[1]): 
           dst[i,j] = src[i,j]

def default_bps(model, c, r, t):
    noBrPointsPerEpoch = model.nbreakpoints
    nleaves = model.nleaves
    nepochs = len(noBrPointsPerEpoch)
    ebps = []
    bps = []
    for e in xrange(0, nepochs):
        theta = 1.0 / c[e]
        nbps = noBrPointsPerEpoch[e]
        if e == nepochs - 1:
            new_bps = [(x*theta+t[e]) for x in expon.ppf([float(i)/nbps for i in xrange(nbps)])]
        else:
            new_bps = linspace(t[e], t[e+1], nbps+1)[:nbps]
        ebps.append(new_bps)
        bps.extend(new_bps)
    return bps, ebps


def logLikelihood(model, obs, col_map, c, r, m, t, posterior_decoding=False):
    noBrPointsPerEpoch = model.nbreakpoints
    nleaves = model.nleaves
    nepochs = len(noBrPointsPerEpoch)
    all_time_breakpoints, time_breakpoints = default_bps(model, c, r, t)

    M = []
    for e in xrange(len(noBrPointsPerEpoch)):
        newM = identity(nleaves)
        newM[:] = m[e]
        M.append(newM)

    pi, T, E = model.run(r, c, time_breakpoints, M, col_map=col_map)
    assert not any(isnan(pi))
    assert not any(isnan(T))
    assert not any(isnan(E))
    logL = inline_forward_scaled(
            array((pi), dtype=float32),
            array((T), dtype=float32),
            array((E), dtype=float32), obs)
    assert logL == logL
    return logL

    # pi_, T_, E_ = model.run(r, c, time_breakpoints, M, col_map=col_map)
    # print "called."
    # assert not any(isnan(pi_))
    # assert not any(isnan(T_))
    # assert not any(isnan(E_))
    # T_ = T_.transpose()
    # E_ = E_.transpose()

    # k = T_.shape[0]
    # L = obs.len()

    # pi = HMMVector(k)
    # T = HMMMatrix(*T_.shape)
    # E = HMMMatrix(*E_.shape)

    # copyTable(T, T_)
    # copyTable(E, E_)
    # for i,v in enumerate(pi_):
    #     pi[i] = v

    # hmm = HMM(pi, T, E)
    # F = HMMMatrix(L,k)
    # scales = HMMVector(L)
    # hmm.forward(obs, scales, F)
    # logL = hmm.likelihood(scales)
    # assert logL == logL
    # if posterior_decoding:
    #     B = HMMMatrix(L,k)
    #     hmm.backward(obs, scales, B)
    #     PD = HMMMatrix(L,k)
    #     hmm.posterior_decoding(obs, F, B, scales, PD)
    #     #pi_count = HMMVector(k)
    #     #T_count = HMMMatrix(*T_.shape)
    #     #E_count = HMMMatrix(*E_.shape)
    #     #hmm.baum_welch(obs, F, B, scales, pi_count, T_count, E_count)
    #     return logL, PD, all_time_breakpoints #, pi_, pi_count, T_count, E_count
    # else:
    #     return logL

# Should not be called directly, create your own version with the params
# needed.
def optimize(file_in, model, seqnames, init, cb=None):
    assert False
    obs, col_map = readObservations(file_in, seqnames)
    nbps = model.nbreakpoints
    logL = [None, _logLikelihood_2, _logLikelihood_3][len(nbps) - 1]
    return fmin(lambda x: -logL(model, obs, *x), init, callback=cb)

