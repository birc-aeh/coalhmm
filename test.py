from scipy import *
from scipy.linalg import expm
from sets import ImmutableSet as iset

#from time_plot import *
from model import Model
from fasta_parser import readAlignment
from hmm import logLikelihood

theta = 2*30000.0 * 25 * 1e-9
#rho = 20000.0 * 0.01 / 1e6

C = 1.0 / theta
R = 1.5e-8 / 1.0e-9

model = Model(2, [5])

#pi, T, E, Q = model.run(rho * (1.0 / theta), 1.0 / theta)
#print E.shape

#time_plot(model.G, Q, 0.004)


alignments = readAlignment("/users/tth/work/simpipe2/3/seq0.fasta")

left, right = alignments["'0'"], alignments["'1'"]
to_char = ['A','C','G','T']
to_val = {'A':0,'C':1,'G':2,'T':3}
left = [to_val[x] for x in left[50000:60000]]
right = [to_val[x] for x in right[50000:60000]]

likelihoods = []
Cs = linspace(0.25*C, 4.0*C)
for c in Cs:
    pi, T, E = model.run(R, c, [[0.0] + [x*theta for x in [.5,1,2,3]]])
    logL = logLikelihood(left, right, T, pi, E)
    print "logL(R=%f, C=%f) = %f" % (R, c, logL)
    likelihoods.append(logL)
from pylab import *
plot(Cs, likelihoods)
vlines(C, min(likelihoods), max(likelihoods))
show()

