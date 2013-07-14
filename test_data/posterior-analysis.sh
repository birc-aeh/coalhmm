#!/bin/sh

for chunk in {0..9};
do
    python ../scripts/posterior-from-I-model.py -T 0.004 -C 497 -R 0.2 analysis/sim_I_${chunk}.ziphmm | \
    python ../scripts/marginal-posterior.py > marginal-posterior-${chunk}.txt

    python ../scripts/posterior-from-I-model.py -T 0.004 -C 497 -R 0.2 --intervals=20 analysis/sim_I_${chunk}.ziphmm | \
    python ../scripts/marginal-posterior.py > marginal-posterior-20states-${chunk}.txt

#    tee posterior-${chunk}.txt | \
done