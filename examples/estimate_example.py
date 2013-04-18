# In this example we will see how to estimate multiple parameters based on an
# alignment.

# We will be using one of two models. The isolation model ('I') or the
# isolation-with-migration model ('IM').
current_model = 'I'

# A model in coalhmm is described by multiple epochs. The first epoch always
# represents the present day and has some number of branches that are completely
# distinct from each other.
# The rest of the epochs can either have a reduced number of branches or keep
# all the previous ones. Additionally it can add migration between branches.
# The reason you might want to keep all branches is to change a parameter over
# time.

# To specify this we need the number of branches in the present and for each
# epoch we need the migrations (if any) and for each epoch transition we must
# know what happens to the branches.

# The transitions are specified as an array where index i tells us which
# branch, branch i is part of in the next epoch.
# If we had three branches and wanted to merge the first and last it would look
# like this: [0,1,0].

# The migration is specified in a similar way, but instead of a single value
# each entry in the array is another array telling us what the branch can
# migrate to.
# If we have three branches and the middle one can migrate to either of its
# neighbours but not vice-versa it would look like this: [[],[0,2],[]].

# Finally coalhmm uses discretized time so we need to specify how finely each
# epoch is subdivided. The present won't benefit from more than 1 time interval
# since everything is separate.
intervals_per_epoch = 10

from coalhmm.model import build_epoch_separated_model
if current_model == 'I':
    modelI = build_epoch_separated_model(
            # With the isolation model we have two branches in the present:
            2,
            # We want to merge 0 and 1, so we map them both to 0:
            [[0,0]],
            [1,intervals_per_epoch])
elif current_model == 'IM':
    # Creating an isolation-with-migration model is similar to the plain
    # isolation model. We add an extra epoch (maintaining the two separate
    # branches) and allow migration from 0 to 1 and from 1 to 0 in it.
    M_for_both = [[1],[0]]
    modelIM = build_epoch_separated_model(
            2,
            [[0,1], [0,0]],
            [1,intervals_per_epoch,intervals_per_epoch],
            # None for no migration in the first and last epoch
            [None, M_for_both, None])

# Now we need some data to work on. 'example_data.fa' is a small chunk of an
# alignment with human, bonobo, chimpanzee and orangutan.
# Lets use the human and chimp.
names = ["hg18", "pantro2"]

# coalhmm comes with a very simple parser for fasta files - you can of course
# also use SeqIO from BioPython or your favorite library.
from coalhmm.fasta_parser import readAlignment

# Once an alignment is loaded we need to convert it to a series of observations
# in the form of a numpy/scipy array of integers.
# In this example we are looking at two species so will use a jukes-cantor like
# model where 0 means no difference, 1 means a difference and 2 means missing
# data.
# We could also have used AA = 0, AC = 1, ..., TT = 15, N* = *N = 16 or any
# other mapping.
from scipy import zeros, int32
def read_observations(filename, names):
    '''Reads in a single fasta file and converts it to a series of
    observations.  The observations assume JC69, so we only care if there is a
    difference or not (or a N).'''
    alignments = readAlignment(filename)
    srcs = [alignments[name] for name in names]
    # alignments = SeqIO.to_dict(SeqIO.parse(open(filename),'fasta'))
    # srcs = [alignments[name].seq for name in names]

    clean = set(['A', 'C', 'G', 'T'])
    A = srcs[0]
    B = srcs[1]
    assert len(A) == len(B)
    L = len(A)
    # int16/int8 can save some memory when using larger datasets
    obs = zeros(L, dtype=int32) 
    for i in xrange(L):
        s1,s2 = A[i], B[i]
        if s1 not in clean or s2 not in clean:
            v = 2
        elif s1 == s2:
            v = 0
        else:
            v = 1
        obs[i] = v
    return obs

all_obs = []
total_L = 0
for seq in ["example_data.fa"]:
    obs = read_observations(seq, names)
    total_L += len(obs)
    all_obs.append(obs)
print "%.2fMbp of data in total" % (total_L/1e6)

# So, now we have our model and our data - how do we estimate something?
# First we need some initial parameters:
g = 20   # generation time
u = 1e-9 # mutation rate
# initial values for recombination rate, population size/coalescence rate and
# migration rate:
i_r = 0.8
i_N = 20e3
i_c = 1.0/(2*g*u*i_N)
i_m = 0.1

# For the I model we guess that the human-chimp split is 5 million years ago.
init_I =   (5.0e6*u, i_c, i_r)
# For the IM model we guess the same, but with an 1 million year of migration.
init_IM =  (4.0e6*u, 5.0e6*u, i_m/2*i_c, i_c, i_r)

# When calculating the emission probabilities coalhmm needs a 'column map' that
# tells it what symbols a value in the observations corresponds to.
# The map is specified in the reverse order though, so it's a mapping from a
# column of symbols to the value that will appear in the observation sequence.
# Since an inverse of the map is made it must be a 1-to-1 mapping.
# In this example we assume JC69 so we have three values, one for two equal
# symbols, one for two different symbols and one for unknown data.
FIXED_COL_MAP = { ('A','A'): 0, ('A','C'): 1, ('N','N'): 2 }

from coalhmm.optimize import logL_multiseq
from scipy.optimize import fmin
def logLikelihood(model, all_obs, c,r,m,t):
    # The simplest way to get a likelihood is to call logL_multiseq which will
    # assign actual times to your intervals, construct the matrices and run the
    # forward algorithm with the resulting matrices.
    return logL_multiseq(model, all_obs, FIXED_COL_MAP, c,r,m,t)

def optimize_f(f_logL, init):
    return fmin(lambda x: -f_logL(*x), init, full_output=True)

# Here are the functions doing the actual estimating of parameters.
# It works by giving scipy.optimize.fmin a function that evaluates the logL
# given the current parameters.
def estimate_I(model, all_obs, T, C, R):
    def logL_all(t, c, r):
        if min([t,c,r]) <= 0:
            return -1e18
        return logLikelihood(model, all_obs, [c]*2, [r]*2, [0]*2, [0.0,t])
    est, L, _, _, _ = optimize_f(logL_all, (T,C,R))
    return (L, list(est))

def estimate_IM(model, all_obs, T1, T2, M, C, R):
    def logL_all(t1, t2, m, c, r):
        if min([t1,t2,m,c,r]) <= 0 or t2 <= t1:
            return -1e18
        return logLikelihood(model, all_obs, [c]*3, [r]*3, [0.0,m,0.0], [0.0,t1,t2])
    est, L, _, _, _ = optimize_f(logL_all, (T1,T2,M,C,R))
    return (L, list(est))

# Now we just need to call one of our functions
if current_model == 'I':
    L, est = estimate_I(modelI, all_obs, *init_I)
elif current_model == 'IM':
    L, est = estimate_IM(modelIM, all_obs, *init_IM)

print 'Final results:'
print "\t".join(map(str, [L] + est))
