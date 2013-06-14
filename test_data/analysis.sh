#!/bin/bash

echo "This script gives an example analysis of both isolation and"
echo "isolation-with-migration models."
echo
echo "It first converts all alignments to the zipHMM format and then"
echo "uses the analysis scripts to estimate parameters."
echo

mkdir analysis

echo "Converting isolation model files..."
I_alignments=`ls I/*.gz`
for algn in $I_alignments; do
    out=`basename $algn .gz`
    cmd="python ../scripts/prepare-alignments.py $algn fasta analysis/${out}.ziphmm"
    echo $cmd
    $cmd
done
echo "Done."
echo


echo "Converting isolation-with-migration model files..."
IM_alignments=`ls IM/*.gz`
for algn in $IM_alignments; do
    out=`basename $algn .gz`
    cmd="python ../scripts/prepare-alignments.py $algn fasta analysis/${out}.ziphmm"
    echo $cmd
    $cmd
done
echo "Done."
echo

echo "Analysis I data with I model."
cmd="python ../scripts/estimate-using-I-model.py --verbose --header --out=I-results-I.txt analysis/sim_I_*"
echo $cmd
$cmd
echo

echo "Analysis I data with IM model."
cmd="python ../scripts/estimate-using-IM-model.py --verbose --header --out=I-results-IM.txt analysis/sim_I_*"
echo $cmd
$cmd
echo
