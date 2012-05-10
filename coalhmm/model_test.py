from model import build_epoch_seperated_model, build_simple_model
from time import time
from scipy import asarray

theta = 2*30000.0 * 25 * 1e-9
C = 1.0 / theta
R = 1.5e-8 / 1.0e-9

def test(nleaves, mappings, nbps, simple=False, migrations=None):
    if simple:
        m = build_simple_model(nleaves, nbps[0])
        print "simple    ",
    else:
        m = build_epoch_seperated_model(nleaves, mappings, nbps, migrations)
        print "seperated ",
    epoch_bps = []
    offset = 0
    for n in nbps:
        epoch_bps.append([(offset + c)*0.5*theta for c in xrange(n)])
        offset += n
    t0 = time()
    pi, T, E = m.run(R, C, epoch_bps)
    pi, T, E = m.run(R, C, epoch_bps)
    pi, T, E = m.run(R, C, epoch_bps)
    t1 = time()
    print "["+str(round(t1-t0, 2)).ljust(5)+"s]",
    if migrations != None:
        print "+flow", str((nleaves, mappings, nbps)).ljust(29), "=>",
    else:
        print str((nleaves, mappings, nbps)).ljust(34), "=>",
    return sum(sum(asarray(E))), sum(pi), sum(sum(T)), T.shape

print test(2, [[0,0]], [1,3])
print test(2, [], [4], True)
print test(3, [[0,0,0]], [1,3], migrations=[[[1],[0],[]], None])
print test(3, [], [4], True)
print test(2, [], [3])
print test(3, [], [3])
print test(2, [[0, 0]], [2, 3])
print test(2, [[0, 1], [0, 0]], [2, 3, 3])
print test(2, [[0, 0], [0, 0]], [2, 3, 3])
print test(2, [], [6], True)
print test(2, [[0, 0], [0, 0]], [2, 5, 3])
print test(2, [[0, 0], [0, 0]], [2, 10, 3])
print test(3, [[0, 0, 0]], [2, 3])
# print test(3, [[0, 0, 1], [0, 0]], [2, 3, 4])

