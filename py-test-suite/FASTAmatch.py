# A chopped up version of repeats.py which counts the number of unique kmers.
# This is done to verify more complicated programs' output.

# For the purposes of this code, a "match" is a specific instance of a type of repeat, while a "repeat" is a type of repeat
# Additionally, "category" is abbreviated to "cat"

from pprint import pprint
from collections import defaultdict
from operator import itemgetter
import itertools
import sys
import os

from testfuncs import *


class Match:

    # An idiosyncrasy that isn't explicitly mentioned in documentation:
    #   For complementary matches, the start position in the repeat is parenthesized.
    #       The start and end positions and the "remaining" value must be changed if the complementary match is to be treated identically to +-strand matches.
    #   For +-stand matches, the "remaining" value in reference to the repeat sequence is parenthesized.
    #   The "remaining" value in reference to the genome is always parenthesized.

    def __init__(self, line, id):
        # self.swScore
        # self.pctdiv
        # self.pctdel
        # self.pctins
        # self.chromName  -  the reference genome file this match came from (typically the chromosome)
        # self.chromStart  -  the starting index in the reference genome
        # self.chromEnd  -  the ending index (exclusive) in the reference genome
        # self.chromRemains
        # self.isComplement  -  the match may be for the complement of the reference sequence
        # self.repeatName  -  The repeat's name (equivalent to its "species" or its leaf in the tree)
        # self.repeatCat  -  the heirarchy of categories this match's repeat belongs to (slash-separated)
        # self.repeatStart  -  the starting index in the repeat consensus sequence
        # self.repeatEnd  -  the ending sequence (exclusive) in the repeat consensus sequence
        # self.repeatRemains
        # self.repeatID - A numerical ID for the repeat type (starts at 1)

        # a unique identifier we add
        self.id = id

        (self.swScore, self.pctdiv, self.pctdel,
         self.pctins, self.chromName, self.chromStart,
         self.chromEnd, self.chromRemains, self.isComplement,
         self.repeatName, self.repeatCat, self.repeatStart,
         self.repeatEnd, self.repeatRemains, self.repeatID) = line.split()

        self.isComplement = self.isComplement == 'C'

        # int-ize the reference coordinates
        self.chromStart, self.chromEnd = int(self.chromStart) - 1, int(self.chromEnd)
        assert self.chromRemains[0] == '(' and self.chromRemains[-1] == ')'
        self.chromRemains = int(self.chromRemains[1:-1])

        # we have to do some converting here to make repeat alignment values in reference to +-strand
        if self.isComplement:
            assert self.repeatStart[0] == '(' and self.repeatStart[-1] == ')'
            self.repeatStart = self.repeatStart[1:-1]
            self.repeatStart, self.repeatEnd, self.repeatRemains = self.repeatRemains, self.repeatStart, self.repeatStart
        else:
            assert self.repeatRemains[0] == '(' and self.repeatRemains[-1] == ')'
            self.repeatRemains = self.repeatRemains[1:-1]

        self.repeatStart, self.repeatEnd = int(self.repeatStart) - 1, int(self.repeatEnd)
        self.repeatRemains, self.repeatID = int(self.repeatRemains), int(self.repeatID)


class Chrom:

    def __init__(self, name):
        self.name = name
        self.seqs = dict()


class Genome:

    # optional kmer library file to prevent need to recompute every time
    def __init__(self, k, genomeName):
        # we begin by loading the full reference genome
        # we store each chromosome as a Chrom instance, which contains a dict of the chroms sequences
        # yes, a defaultdict would be easier, but you can't nest dicts in Python
        self.chroms = dict()
        chromFilenames = filter(lambda s: s.endswith(".fa"), os.listdir(genomeName))
        chromFilepaths = [genomeName + "/" + filename for filename in chromFilenames]
        for filepath in chromFilepaths:
            with open(filepath, "r") as chromFile:
                # take only the filename (e.g. "chr2L.fa") and then drop the ".fa"
                chromName = filepath.split('/')[-1][:-3]
                self.chroms[chromName] = Chrom(chromName)
                for line in chromFile:
                    if line[0] == '>':
                        seqName = line.strip()[1:]
                        # seqs should never be redefined
                        assert seqName not in self.chroms[chromName].seqs
                        # would be more efficient to use bytearrays, but they're unhashable
                        self.chroms[chromName].seqs[seqName] = bytearray()
                    else:
                        self.chroms[chromName].seqs[seqName] += line.strip()

        for key in self.chroms.keys():
            print key

        with open(genomeName + "/" + genomeName + ".fa.out", "r") as matchFile:
            # skip the three header lines
            matchFile.readline()
            matchFile.readline()
            matchFile.readline()
            # generator returning consecutive ints starting at 0
            idGen = itertools.count()
            self.matches = [Match(line, idGen.next()) for line in matchFile]

        # associates kmer with a tuple: (minimizer offset, is reverse complement, LCA repeat ID)
        self.kmers = set()
        numMatches = 0
        print len(self.matches), "matches to process"
        for match in self.matches:
            seq = str(self.chroms[match.chromName].seqs[match.chromName][match.chromStart : match.chromEnd]).lower()
            kmerGen = (min(seq[i:i+k], str(revComp(seq[i:i+k]))) for i in xrange(len(seq) - k + 1))
            for kmer in kmerGen:
                if 'n' not in kmer:
                    self.kmers.add(kmer)
            numMatches += 1
            if numMatches % 10000 == 0:
                print numMatches / 1000, "thousand matches processed"

        print "num unique kmers:", len(self.kmers)


if __name__ == "__main__":
    if len(sys.argv) != 2:
        raise Exception("must supply genome name as command line argument")
    genome = Genome(31, sys.argv[1])
