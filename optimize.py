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

Ne = 25000.0
gen = 20
theta = 2*Ne * gen * 1e-9

C = 1.0 / theta
R = (1.5e-8/gen) / 1.0e-9
tau = 325000e-9

noBrPointsPerEpoch = None
len_obs, obs = 0, None

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
    return len(first), obs

def copyTable(dst, src):
   for i in xrange(src.shape[0]):
       for j in xrange(src.shape[1]): 
           dst[i,j] = src[i,j]

def logLikelihood(model, c, r, t):
    if c < 0 or r < 0 or t < 0:
        return -1e18
    theta = 1.0 / c
    time_breakpoints = [
            [0.0],
            [(x*theta+t) for x in expon.ppf([float(i)/noBrPointsPerEpoch[1] for i in xrange(noBrPointsPerEpoch[1])])]
            ]

    pi_, T_, E_ = model.run(r, c, time_breakpoints)
    T_ = T_.transpose()
    E_ = E_.transpose()
    
    k = T_.shape[0]
    L = len_obs
    
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

def computeLikelihoodProfiles(outputFilePrefix):
   #Ts = linspace(0.125*tau, 10.0*tau)
   #Tlikelihoods = [logLikelihood(C, R, t) for t in Ts]
   Cs = linspace(0.0001*C, 2.0*C)
   Clikelihoods = [logLikelihood(c, R, tau) for c in Cs]

   if outputFilePrefix is not None:
      f = open('%s-T-profile.txt' % outputFilePrefix, 'w')
      print >> f, 'T logL'
      for val,logL in izip(Cs,Clikelihoods):
         print >> f, val, logL
      f.close()

def optimize(file_in, model, seqnames, cb=None):
    global noBrPointsPerEpoch,len_obs, obs
    len_obs, obs = readObservations(file_in, seqnames)
    return fmin(lambda x: -logLikelihood(model, *x), (C,R,tau), callback=cb)

if __name__ == "__main__":
    global model, noBrPointsPerEpoch
    #len_obs, obs = readObservations(sys.argv[1], ["'0'","'1'"])
    #noBrPointsPerEpoch = [1, 3]
    #model = build_epoch_seperated_model(2, [[0,0]], noBrPointsPerEpoch)
    #print optimize(*sys.argv[1:])
    #computeLikelihoodProfiles(sys.argv[2])

    def popsize_to_C(ne):
        return 1.0 / (2.0 * gen * ne * 1e-9)
    def C_to_popsize(c):
        return 1.0 / c / (2 * gen * 1e-9)
    def recrate_to_R(rec):
        return rec / (gen * 0.1)
    def R_to_recrate(r):
        return r * gen * 0.1

    samples = 40
    Cs = map(popsize_to_C, linspace(0.1*Ne, 2.0*Ne, samples))
    Rs = map(recrate_to_R, linspace(0.1, 3.0, samples))
    Ts = linspace(0.5*tau, 5.0*tau, samples)

    #len_obs, obs = readObservations(sys.argv[1], ["'0'","'1'"])
    folder = "self_sim6"
    for nstates in [int(sys.argv[2])]:
        print "Starting the", nstates, "state simulations"
        noBrPointsPerEpoch = [1, nstates]
        model = build_epoch_seperated_model(3, [[0,0,0]], noBrPointsPerEpoch)
        f = open("%s/%s_%i_values.txt" % (folder,sys.argv[1],nstates), 'w')
        print >>f, "C\tR\ttau\tlogL"
        for run in xrange(20):
            filename = "simulated_3s/sim_%s_%i.txt" % (sys.argv[1], run)
            while not os.path.exists(filename):
                time.sleep(10.0)
            c,r,t = optimize(filename, model, ["'0'","'1'","'2'"])
            print >>f, '\t'.join(map(repr, [c,r,t,logLikelihood(model,c,r,t)]))
            f.flush()
        f.close()

        # logFile = open('%s/mle_t_%i_%.txt' % (folder,nstates),'w')
        # print >>logFile, '\t'.join(['Recrate', 'Popsize', 'tau', 'logL'])
        # for t in Ts:
        #     for c in [C]:
        #         for r in [R]:
        #             print >>logFile, '\t'.join(map(repr,
        #                 [R_to_recrate(r), C_to_popsize(c), t, logLikelihood(c, r, t)]))
        # logFile.close()
