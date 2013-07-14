#!/bin/sh

for chunk in {0..9};
do
    python ../scripts/posterior-from-I-model.py -T 0.004 -C 500 -R 0.2 analysis/sim_I_${chunk}.ziphmm | \
    tee posterior-${chunk}.txt | \
    python ../scripts/marginal-posterior.py > marginal-posterior-${chunk}.txt
done