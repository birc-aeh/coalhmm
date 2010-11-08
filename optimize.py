from scipy import *
from scipy.stats import expon
from model import build_simple_model, build_epoch_seperated_model
from fasta_parser import readAlignment
from pyhmmlib import *
from scipy.optimize import fmin
import sys
import os.path
import time
from itertools import izip

def readObservations(filename, seq_names):
    alignments = readAlignment(filename)

    first = alignments[seq_names[0]]
    to_val = {'A':0,'C':1,'G':2,'T':3}
    cols = [[] for _ in range(len(first))]
    for i, seq_name in enumerate(seq_names):
        seq = alignments[seq_name]
        for n in xrange(len(seq)):
            v = seq[n-1]
            cols[n].append(v)
    col_map = dict()
    obs = Sequence(len(first))
    for i, col in enumerate(cols):
        v = col_map.setdefault(tuple(col), len(col_map))
        obs[i] = v
    return obs, col_map

def copyTable(dst, src):
   for i in xrange(src.shape[0]):
       for j in xrange(src.shape[1]): 
           dst[i,j] = src[i,j]

# Needs a few wrappers, since fmin won't work with array inputs
def _logLikelihood_2(model, obs, c, r, m, t):
    if c < 0 or r < 0 or t < 0 or m < 0:
        return -1e18
    return logLikelihood(model, obs, None, [c, c], [r, r], [m, m], [0.0, t])
def _logLikelihood_3(model, obs, c, r, m, t1, t2):
    if c < 0 or r < 0 or t1 < 0 or t2 < t1 or m < 0:
        return -1e18
    return logLikelihood(model, obs, None, [c, c, c], [r, r, r], [m, m, m], [0.0, t1, t2])

def logLikelihood(model, obs, col_map, c, r, m, t, posterior_decoding=False):
    noBrPointsPerEpoch = model.nbreakpoints
    nleaves = model.nleaves
    nepochs = len(noBrPointsPerEpoch)
    time_breakpoints = []
    all_time_breakpoints = []
    for e in xrange(0, nepochs):
        theta = 1.0 / c[e]
        nbps = noBrPointsPerEpoch[e]
        if e == nepochs - 1:
            new_bps = [(x*theta+t[e]) for x in expon.ppf([float(i)/nbps for i in xrange(nbps)])]
        else:
            new_bps = linspace(t[e], t[e+1], nbps+1)[:nbps]
        time_breakpoints.append(new_bps)
        for bp in new_bps:
            all_time_breakpoints.append(bp)

    M = []
    for e in xrange(len(noBrPointsPerEpoch)):
        newM = identity(nleaves)
        newM[:] = m[e]
        M.append(newM)

    pi_, T_, E_ = model.run(r, c, time_breakpoints, M, col_map=col_map)
    assert not any(isnan(pi_))
    assert not any(isnan(T_))
    assert not any(isnan(E_))
    T_ = T_.transpose()
    E_ = E_.transpose()

    k = T_.shape[0]
    L = obs.len()

    pi = HMMVector(k)
    T = HMMMatrix(*T_.shape)
    E = HMMMatrix(*E_.shape)

    copyTable(T, T_)
    copyTable(E, E_)
    for i,v in enumerate(pi_):
        pi[i] = v

    hmm = HMM(pi, T, E)
    F = HMMMatrix(L,k)
    scales = HMMVector(L)
    hmm.forward(obs, scales, F)
    logL = hmm.likelihood(scales)
    assert logL == logL
    if posterior_decoding:
        B = HMMMatrix(L,k)
        hmm.backward(obs, scales, B)
        PD = HMMMatrix(L,k)
        hmm.posterior_decoding(obs, F, B, scales, PD)
        return logL, PD, all_time_breakpoints, pi_
    else:
        return logL

# Should not be called directly, create your own version with the params
# needed.
def optimize(file_in, model, seqnames, init, cb=None):
    assert False
    obs, col_map = readObservations(file_in, seqnames)
    nbps = model.nbreakpoints
    logL = [None, _logLikelihood_2, _logLikelihood_3][len(nbps) - 1]
    return fmin(lambda x: -logL(model, obs, *x), init, callback=cb)

if __name__ == "__main__":
    #len_obs, obs = readObservations(sys.argv[1], ["'0'","'1'"])
    #noBrPointsPerEpoch = [1, 3]
    #model = build_epoch_seperated_model(2, [[0,0]], noBrPointsPerEpoch)
    #print optimize(*sys.argv[1:])
    #computeLikelihoodProfiles(sys.argv[2])

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


    #len_obs, obs = readObservations(sys.argv[1], ["'0'","'1'"])
    folder = "self_sim_3epochs_mig"
    for nstates in [int(sys.argv[2])]:
        print "Starting the", nstates, "state simulations"
        noBrPointsPerEpoch = [1, nstates, nstates]
        model = build_epoch_seperated_model(2, [[0,1],[0,0]], noBrPointsPerEpoch)
        f = open("%s/%s_%i_values_nomig.txt" % (folder,sys.argv[1],nstates), 'w')
        print >>f, "C\tR\tM\ttau_1\ttau_2"
        for run in xrange(20):
            filename = "simulated_3epochs_mig/sim_%s_%i.txt" % (sys.argv[1], run)
            while not os.path.exists(filename):
                time.sleep(10.0)
            c,r,m,t1,t2 = optimize(filename, model, ["'0'","'1'"], (C,R,0.01,0.25*tau,tau))
            print >>f, '\t'.join(map(repr, [c,r,m,t1,t2]))
            print c, r, m, t1, t2
            f.flush()
        f.close()

    #samples = 40
    #Cs = map(popsize_to_C, linspace(0.1*Ne, 2.0*Ne, samples))
    #Rs = map(recrate_to_R, linspace(0.1, 3.0, samples))
    #Ts = linspace(0.5*tau, 5.0*tau, samples)
    #
    # logFile = open('%s/mle_t_%i_%.txt' % (folder,nstates),'w')
    # print >>logFile, '\t'.join(['Recrate', 'Popsize', 'tau', 'logL'])
    # for t in Ts:
    #     for c in [C]:
    #         for r in [R]:
    #             print >>logFile, '\t'.join(map(repr,
    #                 [R_to_recrate(r), C_to_popsize(c), t, logLikelihood(c, r, t)]))
    # logFile.close()
