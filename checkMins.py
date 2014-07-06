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

            for line in f:

                if numLines % 1000000 == 0:
                    print (numLines / 1000000), "million lines processed"
                numLines += 1

                if line[0] == '>':
                    currMin = line[1:-1]
                else:
                    fields = line.split()
                    if len(fields) != 2:
                        raise Exception("ERROR: malformed kmer line")
                    kmer = fields[0]
                    pyMin = mini(kmer, m)
                    if pyMin != currMin:
                        print
                        print line
                        print "kmer:", kmer
                        print "listed min:", currMin
                        print "Python min:", pyMin
                        print
                        sys.exit(0)
