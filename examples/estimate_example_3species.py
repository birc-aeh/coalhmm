import sys
import os

import coalhmm
from coalhmm.model import build_epoch_seperated_model
from coalhmm.optimize import logL_multiseq, readObservations
from scipy import *
from scipy.optimize import fmin

# In this example the col map will be based on the actual observations in the file
COL_MAP = dict()

def logLikelihood(model, all_obs, c,r,m,t):
    return logL_multiseq(model, all_obs, COL_MAP, c,r,m,t)
            #prepare_matrices=zipHMM_prepare_matrices,
            #single_logL=zipHMM_single_logL)

def optimize_f(model, obs, f_logL, init):
    return fmin(lambda x: -f_logL(model, obs, *x), init, full_output=True)

# Here are the functions doing the actual estimating of parameters.
# It works by giving scipy.optimize.fmin a function that evaluates the logL
# given the current parameters. The logging is implemented in an odd way
# that is probably not necessary for your environment.
def estimate_ILS(model, obs, T1, T2, C, R, outfile="/dev/null"):
    def logL_all(model, all_obs, t1, t2, c, r):
        if min([t1,t2,c,r]) <= 0 or t2 <= t1:
            return -1e18
        res = logLikelihood(model, obs, [c]*3, [r]*3, [0.0]*3, [0.0,t1,t2])
        os.system("echo '%s' >> %s" % ('\t'.join(map(str, [t1,t2,c, r, res])), outfile))
        return res
    os.system("echo '%s' > %s" % ('\t'.join(map(str, ["t1","t2", "c", "r", "logL"])), outfile))
    est, L, _, _, _ = optimize_f(model, obs, logL_all, (T1,T2,C,R))
    #       logL.   estimates
    return (L, list(est))


# generation time and mutation rate
g = 20
u = 1e-9

# initial values for recombination rate, population size/coalescence rate and
# migration rate.
i_r = 0.8
i_N = 20e3
i_c = 1.0/(2*g*u*i_N)
i_m = 0.1

estates = 10

name = "ILS.test"

seqs = ["example_data.fa"]

# Creating an isolation model with three species is done by merging species
# pairwise.
modelILS = build_epoch_seperated_model(
        3,
        # we start with a,b,c.
        [
            # merge a and b to '0', c continues as '1'
            [0,0,1],
            # and now '0' and '1' is merged to a new '0'
            [0,0]
        ],
        [1,estates,estates])
nstates = len(modelILS.tree_map)
names = ["bonobo", "pantro2", "hg18"]
init_ILS = (1.0e6*u, 4.5e6*u, i_c, i_r)

all_obs = []
total_L = 0
for seq in seqs:
    # each sequence creates a new column map extended with any new columns that
    # might be seen.
    obs, colmap = readObservations(seq, names, COL_MAP)
    COL_MAP = colmap
    total_L += len(obs)
    all_obs.append(obs)

print "%.2fMbp of data in total" % (total_L/1e6)

L, est = estimate_ILS(modelILS, all_obs, *init_ILS, outfile="/dev/stdout")


print 'Final results:'
print "\t".join(map(str, [L] + est))
