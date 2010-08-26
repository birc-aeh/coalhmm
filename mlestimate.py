from scipy import *
from scipy.stats import expon
from pylab import *
from scipy.optimize import fmin
from model import build_simple_model
from fasta_parser import readAlignment
from pyhmmlib import *


def readData(fname):
   alignments = readAlignment(fname)

   left, right, righter = alignments["'0'"], alignments["'1'"], alignments["'2'"]
   to_char = ['A','C','G','T']
   to_val = {'A':0,'C':1,'G':2,'T':3}
   left = [to_val[x] for x in left]#[100000:200000]]
   right = [to_val[x] for x in right]#[100000:200000]]
   righter = [to_val[x] for x in righter]#[100000:200000]]

   obs = Sequence(len(left))
   for n in xrange(len(left)):
      obs[n] = left[n-1] + (right[n-1] << 2) + (righter[n-1] << 4)

   return obs, len(left)


def copyTable(dst, src):
   for i in xrange(src.shape[0]):
       for j in xrange(src.shape[1]): 
           dst[i,j] = src[i,j]

def logLikelihood(r,c):

    theta = 1.0 / c
    time_breakpoints = [x * theta for x in 
                        expon.ppf([float(i)/noBrPoints for i in xrange(noBrPoints)])]

    pi_, T_, E_ = model.run(r, c, [time_breakpoints])
    T_ = T_.transpose()
    E_ = E_.transpose()
    
    k = T_.shape[0]
    #L = len(obs)

    #return

    pi = HMMVector(k)
    T = HMMTable(*T_.shape)
    E = HMMTable(*E_.shape)

    copyTable(T, T_)
    copyTable(E, E_)
    for i,v in enumerate(pi_):
        pi[i] = v


    hmm = HMM(pi, T, E)
    f = HMMTable(L,k)
    scales = HMMVector(L)
    hmm.forward(obs, scales, f)
    logL = hmm.likelihood(scales)

    return logL

def callBack(x):
   print >> logFile, ' '.join(map(str,x)),
   print >> logFile, logLikelihood(*x)

theta = 2*30000.0 * 25 * 1e-9
C = 1.0 / theta
R = (1.5e-8/25) / 1.0e-9 / 2 # FIXME: not sure about the last / 2 !!!!

obs, L = readData("seq1.fasta")
logFile = open('seq1-5states-2seqs-mle.txt','w')
noSeqs = 3
noIntervals = 10
noBrPoints = noIntervals
model = build_simple_model(noSeqs, noIntervals)

#fmin(lambda x: -logLikelihood(*x), (R,C), callback=callBack)

logLikelihood(R, C)
