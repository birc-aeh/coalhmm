from scipy import *
from scipy.linalg import expm
from scipy.stats import expon
from pylab import *
from model import build_simple_model
from fasta_parser import readAlignment
from pyhmmlib import *
from itertools import izip

alignments = readAlignment("seq0.fasta")

left, right, righter = alignments["'0'"], alignments["'1'"], alignments["'2'"]
to_char = ['A','C','G','T']
to_val = {'A':0,'C':1,'G':2,'T':3}
left = [to_val[x] for x in left]#[100000:200000]]
right = [to_val[x] for x in right]#[100000:200000]]
righter = [to_val[x] for x in righter]#[100000:200000]]


theta = 2*30000.0 * 25 * 1e-9
#rho = 20000.0 * 0.01 / 1e6

C = 1.0 / theta
R = (1.5e-8/25) / 1.0e-9 / 2 # FIXME: not sure about the last / 2 !!!!

noSeqs = 3
obs = Sequence(len(left))
for n in xrange(len(left)):
   obs[n] = left[n-1] + (right[n-1] << 2) + (righter[n-1] << 4)


noIntervals = 2
noBrPoints = noIntervals
model = build_simple_model(noSeqs, noIntervals)

#time_breakpoints = [x * theta for x in [0, 0.5+1e-5, 1]]
#pi, T, E, Q = model.run(C, R, time_breakpoints)


def copyTable(dst, src):
   for i in xrange(src.shape[0]):
       for j in xrange(src.shape[1]): 
           dst[i,j] = src[i,j]

def logLikelihood(r,c):

    theta = 1.0 / c
    #time_breakpoints = [x * theta for x in linspace(0,2,noBrPoints)]
    time_breakpoints = [x * theta for x in 
                        expon.ppf([float(i)/noBrPoints for i in xrange(noBrPoints)])]
    #time_breakpoints = [x * theta for x in linspace(0,4,noBrPoints)]


    pi_, T_, E_ = model.run(r, c, [time_breakpoints])
    T_ = T_.transpose()
    E_ = E_.transpose()
    
    k = T_.shape[0]
    L = len(left)
    
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

    print "logL(R=%f, C=%f) = %f" % (r, c, logL)
    return logL

def computeLikelihoodProfiles(outputFilePrefix):
   Cs = linspace(0.25*C, 4.0*C)
   Rs = linspace(0.25*R, 4.0*R)

   Clikelihoods = [logLikelihood(R,c) for c in Cs]
   Rlikelihoods = [logLikelihood(r,C) for r in Rs]

   if outputFilePrefix is not None:
      f = open('%s-C-profile.txt' % outputFilePrefix, 'w')
      print >> f, 'C logL'
      for c,logL in izip(Cs,Clikelihoods):
         print >> f, c, logL
      f.close()

      f = open('%s-R-profile.txt' % outputFilePrefix, 'w')
      print >> f, 'R logL'
      for r,logL in izip(Rs,Rlikelihoods):
         print >> f, r, logL
      f.close()


   return Cs, Rs, Clikelihoods, Rlikelihoods

def plotProfiles(outputFilePrefix = None):

   Cs, Rs, Clikelihoods, Rlikelihoods = computeLikelihoodProfiles(outputFilePrefix)

   def whichmax(x):
      return x.index(max(x))

   clf()

   subplot(121)
   title('Likelihood of C')
   plot(Cs, Clikelihoods)
   vlines(C, min(Clikelihoods), max(Clikelihoods))
   vlines(Cs[whichmax(Clikelihoods)], min(Clikelihoods), max(Clikelihoods), color='blue')
   xlabel('C')
   ylabel('logL')

   print 'ML estimate for C:', Cs[whichmax(Clikelihoods)]

   subplot(122)
   title('Likelihood of R')
   plot(Rs, Rlikelihoods)
   vlines(R, min(Rlikelihoods), max(Rlikelihoods))
   vlines(Rs[whichmax(Rlikelihoods)], min(Rlikelihoods), max(Rlikelihoods), color='blue')
   xlabel('R')
   ylabel('logL')

   print 'ML estimate for R:', Rs[whichmax(Rlikelihoods)]

   show()


computeLikelihoodProfiles('seq0-2states-equiprob-3seqs')
