# A simple script used to check the correctness of minimizers in a minimizer file.
# Intended for debugging serpentine hand-optimized minimizer algorithms.

import sys
from testfuncs import *

if __name__ == "__main__":

    if sys.argv < 2:
        print "please supply name(s) of minimizer file(s) you wish to check"
        sys.exit(1)

    for filename in sys.argv[1:]:
        f = open(filename, "r")
        m = 15

        for line in f:
            if line[0] == '>':
                currMin = line[1:-1]
            elif line[:2] == "\t\t" or len(line.split()) < 3:
                pass
            else:
                fields = line.split()
                kmer, minOffset, isRevComp = bytearray(fields[0]), int(fields[1]), (fields[2].strip() == '1')
                testMin = revComp(kmer[minOffset:minOffset+m]) if isRevComp else kmer[minOffset:minOffset+m]
                pyMin, pyOffset, pyRC = miniAll(kmer, m)
                if pyMin != testMin or testMin != currMin:
                    print
                    print line
                    print "kmer:", kmer
                    print "minOffset", minOffset
                    print "isRevComp", isRevComp
                    print "listed min:", currMin
                    print "calculated min:", testMin
                    print "Python min:", pyMin, "\t", pyOffset, pyRC
                    print
                    sys.exit(0)
