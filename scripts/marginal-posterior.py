import sys
from numpy import array, mean, zeros

variables = dict()
posterior_means = list()

if len(sys.argv) == 1 or sys.argv[1] == '-':
    input = '/dev/stdin'
else:
    input = sys.argv[1]
    
with open(input) as f:

    no_posteriors = 0
    posterior_sum = None

    for line in f:
        if line.startswith('##'):
            continue
        elif line.startswith('#'):
            var,val = line[1:].split('=')
            variables[var.strip()] = val.strip()
        else:
            post = map(float,line.split('\t')[2:])
            if posterior_sum is None:
                posterior_sum = zeros(len(post))
            assert len(posterior_sum) == len(post), \
                    "%d != %d"%(len(posterior_sum),len(post))
            posterior_sum += post
            no_posteriors += 1


break_times = array(map(float, variables['time_breaks'].split()))
posterior_means = posterior_sum / no_posteriors

print 'break_times\tposteriors'
for i in xrange(len(break_times)):
    print '%f\t%f' % (break_times[i], posterior_means[i])
