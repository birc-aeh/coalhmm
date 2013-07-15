#!/bin/sh

for R in 0.2 0.3 0.4; do
    for chunk in {0..9};
    do
	python ../scripts/posterior-from-I-model.py -T 0.004 -C 500 -R $R analysis/sim_I_${chunk}.ziphmm | \
	    python ../scripts/marginal-posterior.py > marginal-posterior-${chunk}-R-${R}.txt

    done
done
