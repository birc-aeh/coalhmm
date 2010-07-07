from math import log, exp
import scipy
from emission_matrix import cols_to_index
from scipy.stats import expon


def logLikelihood(left, right, T, pi, E):
    """
T is the transition matrix (k, k)
pi is the initial probabilities (k,)
E is the emission matrix (k, 4**nleaves)
    """
    k = T.shape[0]
    L = len(left)
    # hardwire for now...

    F = scipy.zeros((k,L+1))
    F[:,0] = scipy.log(pi)
    T = scipy.log(T)
    E = scipy.log(E)
    
    for n in xrange(1,L+1):
        for i in xrange(k):
            # TODO: Use more columns
            e = E[i, (left[n-1] + (right[n-1] << 2))]
            s = F[0,n-1] + T[0,i]
            for j in xrange(1,k):
                a = F[j,n-1] + T[j,i]
                s += log(1+exp(a-s))
            F[i,n] = e + s

    logLik = F[0,L]
    for i in xrange(1,k):
        logLik += log(1+exp(F[i,L]-logLik))
    
    return logLik

