# A simple script used to check the correctness of minimizers in a minimizer file.
# Intended for debugging serpentine hand-optimized minimizer algorithms.

import sys
from testfuncs import *

if __name__ == "__main__":

    if sys.argv < 2:
        print "please supply name(s) of minimizer file(s) you wish to check"
        sys.exit(1)

    for filename in sys.argv[1:]:
        with open(filename, "r") as f:
            m = 15
            numLines = 0
            mins = set()
            kmers = set()

            for line in f:

                if numLines % 50000 == 0:
                    print (numLines / 1000), "thousand lines processed"
                numLines += 1

                if line[0] == '>':
                    currMin = line[1:-1]
                    if currMin in mins:
                        print "minimizer", currMin, "listed twice"
                        sys.exit(1)
                    mins.add(currMin)
                else:
                    fields = line.split()
                    if len(fields) != 2:
                        raise Exception("ERROR: malformed kmer line")
                    kmer = fields[0]
                    if kmer in kmers:
                        print "kmer", kmer, "listed twice"
                        sys.exit(1)
                    kmers.add(kmer)
                    pyMin = mini(kmer, m)
                    if pyMin != currMin:
                        print
                        print line
                        print "kmer:", kmer
                        print "listed min:", currMin
                        print "Python min:", pyMin
                        print
                        sys.exit(0)

            print "num unique kmers:", len(kmers)
            print "num unique minimizers:", len(mins)
