from scipy import mat, asmatrix, zeros, sum, identity
from scipy.linalg import expm, eig, inv
from math import log, exp

def _generateRateMatrix(states, transitions, rates):
    n_states = len(states)
    Q = asmatrix(zeros((n_states, n_states)))
    for (src, transition, dst) in transitions:
        Q[src,dst] = rates[transition]
    for i in xrange(n_states):
        row = Q[i,:]
        Q[i,i] = -sum(row)
    return Q


def any(it):
    for b in it:
        if b:
            return True
    return False
class CoalTimeComputer(object):

    def __init__(self, rates, states, transitions, epoch_start):
        self.epoch_start = epoch_start
        self.Q = _generateRateMatrix(states, transitions, rates)
        ## a bit of a hack. Exploits that there is a single lineage in E
        ## states and two in B states ... won't generalize, but doesn't have
        ## to for this paper...
        def B(state):
            return len(state) == 2
        def E(state):
            return len(state) == 1

        self.Bstates = [idx for (state,idx) in states.items() if B(state)]
        self.Estates = [idx for (state,idx) in states.items() if E(state)]
        self.state_range = range(len(states))

        ## Eigenvalue decomposition ... I checked that it is roughly correct
        ## (with some nummerical errors): Q = Zmatrix * diag(psi) * Zm1
        self.psi, self.Zmatrix = eig(self.Q)
        self.Zmatrix = mat(self.Zmatrix)
        self.Zm1 = self.Zmatrix.I


    ## FIXME: Tabulate results for dynamic programming
    def Z(self, a,x,y,b,alpha,beta):
        res = self.Q[x,y] * self.Zmatrix[a,alpha] * self.Zm1[alpha,x] * self.Zmatrix[y,beta] * self.Zm1[beta,b]
        #if res < 0:
        #    print "Z(...) < 0"
        #    print "Q[x,y]       =", self.Q[x,y]
        #    print "Z[a,alpha]   =", self.Zmatrix[a,alpha]
        #    print "Zm1[alpha,x] =", self.Zm1[alpha,x]
        #    print "Z[y,beta]    =", self.Zmatrix[y,beta]
        #    print "Zm1[beta,b]  =", self.Zm1[beta,b]
        #    assert False
        return res

    def I(self, lam, Del):
        if lam == 0.0:
            return Del**2 / 2
        else:
            #return (Del - 1/lam) * exp(Del * lam) / lam + 1/lam**2
            part1 = (Del - 1/lam) * exp(Del * lam) / lam
            part2 = 1/lam**2
            if part1 < 0 and part2 > 0:
                if log(-part1) - log(part2) == 0.0:
                    return 0.0
            if part2 < 0 and part1 > 0:
                if log(-part2) - log(part1) == 0.0:
                    return 0.0
            return part1 + part2

    ## FIXME: Tabulate results for dynamic programming
    def J(self, a,x,y,b,DeltaI):
        return sum(self.Z(a,x,y,b,alpha,beta)*self.I(self.psi[alpha]-self.psi[beta],DeltaI)*exp(self.psi[beta]*DeltaI)
                   for alpha in self.state_range
                   for beta  in self.state_range)

    ## FIXME: Assume that DeltaI is the same for all intervals and
    ## tabulate the M(a) results.
    def M(self, a, DeltaI, PDeltaI):
        assert DeltaI > 0
        # The math says to multiply the inner sum by PDeltaI[a,b]
        return sum(sum(self.J(a,x,y,b,DeltaI) 
                                      for x in self.Bstates
                                      for y in self.Estates)
                   for b in self.Estates)


    def __call__(self, tauI, tauIp1):
        tauI = tauI - self.epoch_start
        tauIp1 = tauIp1 - self.epoch_start
        DeltaI = tauIp1 - tauI

        P0i = expm(self.Q*tauI)
        #if abs(tauI) < 1e-8:
        #    P0i = identity(P0i.shape[0])
        PDeltaI = expm(self.Q*DeltaI)

        PLinI = sum(P0i[0,a]*PDeltaI[a,b]
                for a in self.Bstates
                for b in self.Estates)

        Msum = sum(P0i[0,a]*self.M(a,DeltaI,PDeltaI)
                                    for a in self.Bstates)
        res = tauI + 1.0/PLinI * Msum
        res = res.real
        if res > tauIp1 or res < tauI:
            print "Nooo!", tauI, tauIp1, res, self.epoch_start
            print "  ", PLinI, 1.0/PLinI, Msum.real
            print "  ", 1.0/PLinI * Msum.real
        return self.epoch_start + res
