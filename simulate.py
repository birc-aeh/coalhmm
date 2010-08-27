from scipy import *
from scipy.linalg import expm
from scipy.stats import expon
from model import build_simple_model, build_epoch_seperated_model
from fasta_parser import readAlignment
from pyhmmlib import *
from itertools import izip
from scipy.optimize import fmin
from numpy.random import multinomial
import os
import sys

do_progress_bar = False
if len(sys.argv) > 2:
    do_progress_bar = bool(sys.argv[2])

Ne = 25000.0
gen = 20
theta = 2*Ne * gen * 1e-9

C = 1.0 / theta
R = (1.5e-8/gen) / 1.0e-9
tau = 325000e-9

noBrPointsPerEpoch = [1, 6]
model = None

def choose_weighted(P):
    return argmax(multinomial(1, P))

def index_to_cols(value, n):
    """Converts an integer to a column of n symbols"""
    for i in xrange(n):
        yield i, (value >> 2*i) & 3

def simulate(c,r,t,N):
    theta = 1.0 / c
    time_breakpoints = [
            [0.0],
            [(x*theta+t) for x in expon.ppf([float(i)/noBrPointsPerEpoch[1] for i in xrange(noBrPointsPerEpoch[1])])]
            ]

    print "  Matrices generated"
    pi, T, E = model.run(r, c, time_breakpoints)
    species = model.nleaves
    S = zeros(N, int)
    columns = [zeros(N, int) for _ in xrange(species)]
    E = asarray(E)

    S[0] = choose_weighted(pi)
    if do_progress_bar:
        print "  [" + ' '*100 + ']',
    for i in xrange(1, N):
        progress = 1.0*i/N
        if do_progress_bar:
            print "\r  [" + '='*int(progress*100) + ' '*int((1.0 - progress)*100) + ']',
        col = choose_weighted(E[S[i-1],:])
        for ci, cv in index_to_cols(col, species):
            columns[ci][i-1] = cv
        S[i] = choose_weighted(T[:,S[i-1]])

    if do_progress_bar:
        print "\r" + ' '*104, "\r",

    return columns

folder = "simulated_3s"
char_map = ['A', 'C', 'G', 'T']
for nstates in [int(sys.argv[1])]:
    for run in xrange(20):
        print "Simulating with", nstates, "states"
        print "  C =", C, "R =", R, "tau = ", tau
        noBrPointsPerEpoch = [1, nstates]
        model = build_epoch_seperated_model(3, [[0,0,0]], noBrPointsPerEpoch)
        print "  Model ready"
        cols = simulate(C, R, tau, 1000000)
        print "  Simulation done, writing result"
        filename = '%s/sim_%i_%i.txt' % (folder,nstates,run)
        fastaFile = open(filename+".tmp",'w')
        for i, row in enumerate(cols):
            print >>fastaFile, ">'%i'" % i
            line = ("".join([char_map[v] for v in row]))
            print >>fastaFile, line
        fastaFile.close()
        os.system("mv " + filename+".tmp " + filename)

