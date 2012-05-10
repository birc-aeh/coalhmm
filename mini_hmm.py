from scipy import *
import scipy.weave as weave
import numpy as np

def logaddexp(a, b):
    return a + log(1 + exp(b-a))
X = logaddexp

def inline_forward(Pi, T, E, obs):
    k = len(T)
    L = len(obs)
    A = zeros((L,k), dtype='float64')
    Ew = E.shape[1]

    code = """
    #line 16 "mini_hmm.py"
    #define logaddexp(a,b) ((a) + log(1.0 + exp((b)-(a))))
    int t,j,i;
    double x;
    for (i = 0; i < k; i++)
        A[0*k + i] = Pi[i] + E[obs[0] + i*Ew];
    for (t = 1; t < L; t++)
        for (j = 0; j < k; j++)
        {
            x = T[0*k + j] + A[(t-1)*k + 0];
            for (i = 1; i < k; i++)
                x = logaddexp(x, T[i*k + j] + A[(t-1)*k + i]);
            A[t*k + j] = E[obs[t-1] + j*Ew] + x;
        }
    x = A[(L-1)*k + 0];
    for (i = 1; i < k; i++)
        x = logaddexp(x, A[(L-1)*k + i]);
    return_val = x;
    """

    res = weave.inline(code,
            ['k', 'L', 'A', 'Pi', 'T', 'E', 'obs', 'Ew'],
            compiler="gcc")
    return res

def inline_light_forward(Pi, T, E, obs):
    k = len(T)
    L = len(obs)
    A = zeros((2,k), dtype='float64')
    Ew = E.shape[1]

    code = """
    #line 48 "mini_hmm.py"
    #define logaddexp(a,b) ((a) + log(1.0 + exp((b)-(a))))
    int t,j,i,o;
    double x;
    for (i = 0; i < k; i++)
        A[0*k + i] = Pi[i] + E[obs[0]*Ew + i];
    for (t = 1; t < L; t++)
    {
        o = obs[t-1];
        for (j = 0; j < k; j++)
        {
            x = T[0 + j*k] + A[(t-1 & 1)*k + 0];
            for (i = 1; i < k; i++)
                x = logaddexp(x, T[i + j*k] + A[(t-1 & 1)*k + i]);
            A[(t&1)*k + j] = E[o*Ew + j] + x;
        }
    }
    x = A[(L-1 & 1)*k + 0];
    for (i = 1; i < k; i++)
        x = logaddexp(x, A[(L-1 & 1)*k + i]);
    return_val = x;
    """

    res = weave.inline(code,
            ['k', 'L', 'A', 'Pi', 'T', 'E', 'obs', 'Ew'],
            compiler="gcc")
    return res

