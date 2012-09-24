from scipy import *
import scipy.weave as weave
import numpy as np

def inline_forward_scaled(Pi, T, E, obs):
    k = len(T)
    L = len(obs)
    An = zeros((L,k), dtype='float64')
    C = zeros(L, dtype='float64')
    D = zeros(k, dtype='float64')
    Ew = E.shape[1]

    code = """
    #line 18 "mini_hmm.py"
    int t,j,i,o;
    double x, C_n;
    C_n = 0.0;
    for (i = 0; i < k; i++)
        C_n += Pi[i] * E[obs[0] + i*Ew];
    C[0] = C_n;
    for (i = 0; i < k; i++)
        An[0*k + i] = Pi[i] * E[obs[0] + i*Ew] / C_n;
    for (t = 1; t < L; t++)
    {
        o = obs[t-1];
        for (j = 0; j < k; j++)
        {
            x = 0;
            for (i = 0; i < k; i++)
                x += T[i*k + j] * An[(t-1)*k + i];
            D[j] = x * E[o + j*Ew];
        }
        C_n = 0.0;
        for (j = 0; j < k; j++)
            C_n += D[j];
        C[t] = C_n;
        for (j = 0; j < k; j++)
            An[t*k + j] = D[j]/C_n;
    }
    x = 0.0;
    for (t = 0; t < L; t++)
        x += log(C[t]);
    return_val = x;
    """

    res = weave.inline(code,
            ['k', 'L', 'An', 'C', 'D', 'Pi', 'T', 'E', 'obs', 'Ew'],
            compiler="gcc")
    return res

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
        A[0*k + i] = Pi[i] + E[obs[0] + i*Ew];
    for (t = 1; t < L; t++)
    {
        o = obs[t-1];
        for (j = 0; j < k; j++)
        {
            x = T[0 + j*k] + A[((t-1) & 1)*k + 0];
            for (i = 1; i < k; i++)
                x = logaddexp(x, T[i + j*k] + A[((t-1) & 1)*k + i]);
            A[(t&1)*k + j] = E[o + j*Ew] + x;
        }
    }
    x = A[((L-1) & 1)*k + 0];
    for (i = 1; i < k; i++)
        x = logaddexp(x, A[((L-1) & 1)*k + i]);
    return_val = x;
    """

    res = weave.inline(code,
            ['k', 'L', 'A', 'Pi', 'T', 'E', 'obs', 'Ew'],
            compiler="gcc")
    return res

