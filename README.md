CoalHMM
=======

This is CoalHMM, a framework for demographic inference using a sequential Markov
coalescence model. The framework is based on continuous time Markov chains of systems
of two neighbouring nucleotides and allows for construction of complex demographic
models, only limited by state space explosion (which actually is a rather hard limit
when it comes to the number of populations and samples it can handle).

The framework has been used to construct two concrete models, an isolation model
described in

    T. Mailund, J. Y. Dutheil, A. Hobolth, G. Lunter, and M. H. Schierup, “Estimating
    divergence time and ancestral effective population size of Bornean and Sumatran
    orangutan subspecies using a coalescent hidden Markov model.,” PLoS Genet, vol. 7,
    no. 3, p. e1001319, Mar. 2011.

and an isolation-with-migration model described in

    T. Mailund, A. E. Halager, M. Westergaard, J. Y. Dutheil, K. Munch, L. N. Andersen,
    G. Lunter, K. Prüfer, A. Scally, A. Hobolth, and M. H. Schierup, “A New Isolation
    with Migration Model along Complete Genomes Infers Very Different Divergence
    Processes among Closely Related Great Ape Species.,” PLoS Genet, vol. 8, no. 12,
    pp. e1003125–e1003125, Nov. 2012.

These two models are available in scripts/estimate-using-I-model.py and scripts/estimate-using-IM-model.py, respectively.

Usage
-----

Both models work on pairwise alignments.  Given a pairwise alignment, it must first be
translated into the compressed format used by zipHMM (http://ziphmm.googlecode.com)
that we use for fast likelihood computations.  The script prepare-alignments.py reads
a number of different alignments formats and outputs an alignment in zipHMM format.

If you have full genome alignments you might want to split it into segments that you can
handle in parallel, but the zipHMM should be fast enough to handle full chromosomes if
you prefer.

Once you have formatted the alignment, you can run isolation or isolation-with-migration model using the estimate-using-I-model.py or estimate-using-IM-model.py scripts.
Both scripts accept initial parameters for the likelihood maximisation and will output
the estimated parameters.

See the script test_data/analysis.sh for an example of preparing alignments and running
the two models.

