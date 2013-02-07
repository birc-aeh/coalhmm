from scipy import zeros
import scipy.weave as weave
import numpy as np

def inline_forward_scaled(pi, T, E, obs):
    forward, logL = calc_forward_and_logL(pi, T, E, obs)
    return logL


def calc_forward_and_logL(pi, T, E, obs):
    k = len(T)
    L = len(obs)
    An = zeros((L,k), dtype=np.float64)
    C = zeros(L, dtype=np.float64)
    D = zeros(k, dtype=np.float64)
    Ew = E.shape[1]

    code = """
    #line 18 "mini_hmm.py"
    int t,j,i,o;
    double x, C_n;
    C_n = 0.0;
    for (i = 0; i < k; i++)
        C_n += pi[i] * E[obs[0] + i*Ew];
    C[0] = C_n;
    for (i = 0; i < k; i++)
        An[0*k + i] = pi[i] * E[obs[0] + i*Ew] / C_n;
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
            ['k', 'L', 'An', 'C', 'D', 'pi', 'T', 'E', 'obs', 'Ew'],
            compiler="gcc")
    return An, res

