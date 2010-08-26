from scipy import *
from scipy.linalg import expm
from scipy.stats import expon
from pylab import *
from model import build_simple_model, build_epoch_seperated_model
from fasta_parser import readAlignment
from pyhmmlib import *
from itertools import izip
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

def read_observations(filename):
    alignments = readAlignment(filename)

    left, right, righter = alignments["'0'"], alignments["'1'"], alignments["'2'"]
    to_char = ['A','C','G','T']
    to_val = {'A':0,'C':1,'G':2,'T':3}
    left = [to_val[x] for x in left]
    right = [to_val[x] for x in right]
    righter = [to_val[x] for x in righter]
    obs = Sequence(len(left))
    for n in xrange(len(left)):
       obs[n] = right[n-1] + (righter[n-1] << 2)# + (righter[n-1] << 4)
    return len(left), obs

len_obs, obs = read_observations("../simulated_data/Mig750_1/Sim10/seq.unif.fasta")

noBrPointsPerEpoch = [1, 6]
model = build_epoch_seperated_model(2, [[0,0]], noBrPointsPerEpoch)


def copyTable(dst, src):
   for i in xrange(src.shape[0]):
       for j in xrange(src.shape[1]): 
           dst[i,j] = src[i,j]

def logLikelihood(c,r,t):
    theta = 1.0 / c
    time_breakpoints = [
            [0.0],
            [(x*theta+t) for x in expon.ppf([float(i)/noBrPointsPerEpoch[1] for i in xrange(noBrPointsPerEpoch[1])])]
            #[(x*theta+t) for x in expon.ppf(linspace(0.0, 0.999, noBrPointsPerEpoch[1]))]
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

    print "logL(R=%g, C=%g, t=%g) = %g" % (r, c, t, logL)
    return logL

def computeLikelihoodProfiles(outputFilePrefix):
   #Ts = linspace(0.125*tau, 10.0*tau)
   #Tlikelihoods = [logLikelihood(C, R, t) for t in Ts]
   Cs = linspace(0.0001*C, 2.0*C)
   Clikelihoods = [logLikelihood(c, R, tau) for c in Cs]

   if outputFilePrefix is not None:
      f = open('%s-T-profile.txt' % outputFilePrefix, 'w')
      print >> f, 'T logL'
      for t,logL in izip(Cs,Clikelihoods):
         print >> f, t, logL
      f.close()

   return Ts, Tlikelihoods

samples = 40
Cs = map(popsize_to_C, linspace(0.1*Ne, 2.0*Ne, samples))
Rs = map(recrate_to_R, linspace(0.1, 3.0, samples))
Ts = linspace(0.5*tau, 6.0*tau ,samples)

folder = "2_locked"
for nstates in [2, 3, 4, 5, 6, 7, 8, 9, 10]:
    noBrPointsPerEpoch = [1, nstates]
    model = build_epoch_seperated_model(2, [[0,0]], noBrPointsPerEpoch)
    #model = build_simple_model(2, nstates)
    logFile = open('%s/mle_C_%i.txt' % (folder,nstates),'w')
    print >>logFile, '\t'.join(['Recrate', 'Popsize', 'tau', 'logL'])
    for t in [tau]:
        for c in Cs:
            for r in [R]:
                print >>logFile, '\t'.join(map(repr,
                    [R_to_recrate(r), C_to_popsize(c), t, logLikelihood(c, r, t)]))
    logFile.close()


# for c in [750, 250, 500, 0]:
#     for i in xrange(1,26):
#         print "======="
#         print ("mig%i_1, sim%i" % (c,i))
#         len_obs, obs = read_observations("../simulated_data/Mig%i_1/Sim%i/seq.unif.fasta" % (c,i))
#         logFile = open('mig%i_1-sim%i-2states-mle.txt' % (c,i),'w')
#         def callBack(x):
#             print >> logFile, ' '.join(map(str,x)),
#             print >> logFile, logLikelihood(*x)
#         #computeLikelihoodProfiles('mig750_1-sim1-2states-2seqs')
#         fmin(lambda x: -logLikelihood(*x), (C,R,tau), callback=callBack)
#         logFile.close()
