from scipy import *
from scipy.linalg import expm
from sets import ImmutableSet as iset

from time_plot import *
from model import Model



theta = 20000.0 * 20 * 1e-9
rho = 20000.0 * 0.01 / 1e6

model = Model(2, 5)

pi, T, E, Q = model.run(rho * (1.0 / theta), 1.0 / theta)
print E.shape
print E
#print "lala"
#E, T, pi = model.run(rho * (2.0 / theta), 2.0 / theta)

#time_plot(model.G, Q, 0.004)

