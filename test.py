from scipy import *
from scipy.linalg import expm
from sets import ImmutableSet as iset

#from time_plot import *
from model import build_epoch_seperated_model, build_simple_model
from fasta_parser import readAlignment
from hmm import logLikelihood

theta = 2*30000.0 * 25 * 1e-9
#rho = 20000.0 * 0.01 / 1e6

C = 1.0 / theta
R = 1.5e-8 / 1.0e-9

#nbps = [2, 3, 3]
#model = build_epoch_seperated_model(3, [[0, 0, 1], [0, 0]], nbps)
nbps = [6]
model = build_simple_model(3, nbps[0])

#pi, T, E, Q = model.run(rho * (1.0 / theta), 1.0 / theta)
#print E.shape

#time_plot(model.G, Q, 0.004)


alignments = readAlignment("/users/tth/work/simpipe2/3/seq5.fasta")

left, right = alignments["'0'"], alignments["'1'"]
to_char = ['A','C','G','T']
to_val = {'A':0,'C':1,'G':2,'T':3}
left = [to_val[x] for x in left[5000:60000]]
right = [to_val[x] for x in right[5000:60000]]

likelihoods = []
Cs = linspace(0.25*C, 4.0*C, 10)
for cc in Cs:
    theta = 1/cc
    epoch_bps = []
    offset = 0
    for n in nbps:
        epoch_bps.append([(offset + c)*0.5*theta for c in xrange(n)])
        offset += n
    pi, T, E = model.run(R, cc, epoch_bps)
    logL = logLikelihood(left, right, T, pi, E)
    print "logL(R=%f, C=%f) = %f" % (R, cc, logL)
    likelihoods.append(logL)
from pylab import *
plot(Cs, likelihoods)
vlines(C, min(likelihoods), max(likelihoods))
show()

