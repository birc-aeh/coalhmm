import os, os.path, sys
import gzip
import tempfile
from optparse import OptionParser

from pyZipHMM import Forwarder
from Bio import SeqIO

def main():
    usage="""%prog [options] <input> <input format> <output dir>

This program reads in an input sequence in any format supported by BioPython
and writes out a preprocessed file ready for use with zipHMM.
Also supports gzipped input files, if the name ends with `.gz`.

Assumption #1: Either the file is a pairwise alignment, or you have provided
exactly two names to the `--names` option.

Assumption #2: The file uses a simple ACGT format (and N/-). Anything else will
be interpreted as N and a warning will be given with all unknown symbols.

Warning: This program uses SeqIO.to_dict to read in the entire alignment, you
may want to split the alignment first if it's very large.
"""


    parser = OptionParser(usage=usage, version="%prog 1.0")

    parser.add_option("--names",
                      dest="names",
                      type="string",
                      default=None,
                      help="A comma-separated list of names to use from the source file")
    parser.add_option("-v", "--verbose",
                      dest="verbose",
                      action="store_true",
                      default=False,
                      help="Print some stuff")

    (options, args) = parser.parse_args()

    if len(args) != 3:
        parser.error("Needs input file, input format and output file")
    in_filename = args.pop(0)
    in_format = args.pop(0)
    output_dirname = args.pop(0) 

    assert os.path.exists(in_filename), "Must use an existing input file"
    if in_filename.endswith('.gz'):
        if options.verbose:
            print "Assuming '%s' is a gzipped file." % in_filename
        inf = gzip.open(in_filename)
    else:
        inf = open(in_filename)
    
    if options.verbose:
        print "Loading data...",
        sys.stdout.flush()
    alignments = SeqIO.to_dict(SeqIO.parse(inf,in_format))
    if options.verbose:
        print "done"

    if options.names:
        names = options.names.split(',')
    else:
        names = list(alignments.keys())
    assert len(names) == 2, "Must be a pairwise alignment."
    if options.verbose:
        print "Assuming pairwise alignment between '%s' and '%s'" % (names[0],names[1])
    srcs = [alignments[name].seq for name in names]

    clean = set('ACGT')
    A = srcs[0]
    B = srcs[1]
    assert len(A) == len(B)
    L = len(A)
    fd, foutname = tempfile.mkstemp()
    if options.verbose:
        print "Writing temp file readable by zipHMM to '%s'..." % (foutname),
        sys.stdout.flush()
    seen = set()
    with os.fdopen(fd, 'w', 64*1024) as f:
        for i in xrange(L):
            s1,s2 = A[i].upper(), B[i].upper()
            seen.add(s1)
            seen.add(s2)
            if s1 not in clean or s2 not in clean:
                print >>f, 2,
            elif s1 == s2:
                print >>f, 0,
            else:
                print >>f, 1,
    if options.verbose:
        print "done"
    if len(seen - set('ACGTN-')) > 1:
        print >>sys.stderr, "I didn't understand the following symbols form the input sequence: %s" % (''.join(list(seen - set('ACGTN-'))))
    if options.verbose:
        print "zipHMM  is preprocessing...",
        sys.stdout.flush()
    f = Forwarder.fromSequence(seqFilename = foutname,
                               alphabetSize = 3, minNoEvals = 500)
    if options.verbose:
        print "done"

    if options.verbose:
        print "Writing zipHMM data to '%s'..." % (output_dirname),
        sys.stdout.flush()
    if not os.path.exists(output_dirname):
        os.makedirs(output_dirname)
    f.writeToDirectory(output_dirname)
    os.rename(foutname, os.path.join(output_dirname, 'original_sequence'))
    if options.verbose:
        print "done"

if __name__ == "__main__":
    main()
