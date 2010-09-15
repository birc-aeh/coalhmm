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
    obs = Sequence(len(first))
    for i, seq_name in enumerate(seq_names):
        seq = alignments[seq_name]
        for n in xrange(len(seq)):
            v = to_val[seq[n-1]]
            obs[n] += v << 2*i
    return obs

def copyTable(dst, src):
   for i in xrange(src.shape[0]):
       for j in xrange(src.shape[1]): 
           dst[i,j] = src[i,j]

# Needs a few wrappers, since fmin won't work with array inputs
def _logLikelihood_2(model, obs, c, r, m, t):
    if c < 0 or r < 0 or t < 0 or m < 0:
        return -1e18
    return logLikelihood(model, obs, [c, c], [r, r], [m, m], [0.0, t])
def _logLikelihood_3(model, obs, c, r, m, t1, t2):
    if c < 0 or r < 0 or t1 < 0 or t2 < t1 or m < 0:
        return -1e18
    return logLikelihood(model, obs, [c, c, c], [r, r, r], [m, m, m], [0.0, t1, t2])

def logLikelihood(model, obs, c, r, m, t):
    noBrPointsPerEpoch = model.nbreakpoints
    assert noBrPointsPerEpoch[0] == 1
    nleaves = model.nleaves
    time_breakpoints = [[0.0]]
    for e in xrange(1, len(noBrPointsPerEpoch)):
        theta = 1.0 / c[e]
        nbps = noBrPointsPerEpoch[e]
        time_breakpoints.append(
                [(x*theta+t[e]) for x in expon.ppf([float(i)/nbps for i in xrange(nbps)])]
                )

    M = []
    for e in xrange(len(noBrPointsPerEpoch)):
        newM = identity(nleaves)
        newM[:] = m[e]
        M.append(newM)

    pi_, T_, E_ = model.run(r, c, time_breakpoints, M)
    #print pi_
    #print T_
    #print E_
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
    f = HMMMatrix(L,k)
    scales = HMMVector(L)
    hmm.forward(obs, scales, f)
    logL = hmm.likelihood(scales)
    return logL

def optimize2(file_in, model, seqnames, init, cb=None):
    obs = readObservations(file_in, seqnames)
    nbps = model.nbreakpoints
    logL = [None, _logLikelihood_2, _logLikelihood_3][len(nbps) - 1]
    return fmin(lambda x: -logL(model, obs, x[0], x[1], x[2], 0.25*x[3], x[3]), init, callback=cb)

def optimize(file_in, model, seqnames, init, cb=None):
    obs = readObservations(file_in, seqnames)
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
