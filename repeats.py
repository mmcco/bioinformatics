# complement support needs to be added

# For the purposes of this code, a "match" is a specific instance of a type of repeat, while a "repeat" is a type of repeat
# Additionally, "category" is abbreviated to "cat"

from pprint import pprint
from collections import defaultdict
from operator import itemgetter
import itertools


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
        self.chromStart, self.chromEnd = int(self.chromStart), int(self.chromEnd)
        assert self.chromRemains[0] == '(' and self.chromRemains[-1] == ')'
        self.chromRemains = int(self.chromRemains[1:-1])

        # we have to do some converting here to make repeat alignment values in reference to +-strand
        if self.isComplement:
            assert self.repeatStart[0] == '(' and self.repeatStart[-1] == ')'
            self.repeatStart = self.repeatStart[1:-1]
            self.repeatStart, self.repeatEnd, self.repeatRemains = self.repeatRemains, self.repeatEnd, self.repeatStart
        else:
            assert self.repeatRemains[0] == '(' and self.repeatRemains[-1] == ')'
            self.repeatRemains = self.repeatRemains[1:-1]

        self.repeatStart, self.repeatEnd = int(self.repeatStart), int(self.repeatEnd)
        self.repeatRemains, self.repeatID = int(self.repeatRemains), int(self.repeatID)

        # populated by self.minimize
        self.kmers = []


class Chrom:

    def __init__(self, name):
        self.name = name
        self.seqs = dict()


class Repeat:

    def __init__(self, cat, id):
        self.id = id
        # construct a full ancestry path for this repeat
        self.cat = cat
        # we also store a concatenated version of self.cat
        # a bit of a hack, but constantly joining and splitting takes time, and there aren't many cats (or repeats)
        self.fullName = '/'.join(cat)


# this class represents a node in tree of repeat categories
# it can represent either a repeat or a repeat category
# i.e. for the repeat DNA/P/PROTOP we add CatNodes for DNA, DNA/P, and DNA/P/PROTOP
class CatNode:
    def __init__(self, cat):
        self.cat = cat
        self.fullName = '/'.join(cat)

        if self.fullName == "root":
            self.parentName = None
        elif cat[0] == "Other" or len(cat) == 1:
            self.parentName = "root"
        else:
            self.parentName = '/'.join(cat[:-1])

        # in the future, we may want to store a pointer to the associated Repeat, for CatNodes that represent a Repeat

        # makeTree populates self.parent and self.children
        self.parent = None
        self.children = set()


class Genome:

    # optional kmer library file to prevent need to recompute every time
    def __init__(self, k, m, chromFilepaths, matchFilepath, kmerLibrary=None):
        # we begin by loading the full reference genome
        # we store each chromosome as a Chrom instance, which contains a dict of the chroms sequences
        # yes, a defaultdict would be easier, but you can't nest dicts in Python
        self.chroms = dict()
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
                        self.chroms[chromName].seqs[seqName] = bytearray()
                    else:
                        self.chroms[chromName].seqs[seqName] += line.strip()

        for key in self.chroms.keys():
            print key

        with open(matchFilepath, "r") as matchFile:
            # skip the three header lines
            matchFile.readline()
            matchFile.readline()
            matchFile.readline()
            # generator returning consecutive ints starting at 0
            idGen = itertools.count()
            self.matches = [Match(line, idGen.next()) for line in matchFile]

        # a list of all unique repeat types, indexed by ID
        self.repeats = []
        # maps repeat names to their corresponding repeats
        self.repeatsDict = {}
        for match in self.matches:
            # we only need to add each repeat once
            if len(self.repeats) > match.repeatID and self.repeats[match.repeatID]:
                continue
            else:
                cat = match.repeatCat.split('/')
                cat.append(match.repeatName)
                repeat = Repeat(cat, match.repeatID)
                # we're using the slightly hacky method of growing the list as needed
                while len(self.repeats) < repeat.id + 1:
                    self.repeats.append(None)
                self.repeats[repeat.id] = repeat
                self.repeatsDict[repeat.fullName] = repeat

        # a dict of CatNodes indexed by name
        self.catNodes = {"root": CatNode(["root"])}
        self.tree = self.catNodes["root"]
        for repeat in self.repeatsDict.values():
            # executes for every heirarchical level (see above DNA/P/PROTOP example)
            branchGen = (repeat.cat[:i+1] for i in xrange(len(repeat.cat)))
            for cat in branchGen:
                if '/'.join(cat) not in self.catNodes:
                    catNode = CatNode(cat)
                    self.catNodes[catNode.fullName] = catNode

        # we now loop through all CatNodes to populate the parents and children attrs
        for catNode in self.catNodes.values():
            if catNode.parentName:
                catNode.parent = self.catNodes[catNode.parentName]
                catNode.parent.children.add(catNode)

        # length of each k-mer
        self.k = k
        # length of each minimizer
        self.m = m
        # populated later by minimize()
        self.kmers = defaultdict(list)

        # kmers stored in a tuple: (kmer, minimizer offset, chromName)
        if kmerLibrary:
            self.parseKmerLibrary(kmerLibrary)
        else:
            self.minimize()
        # use Kraken technique, sorting primarily on minimizer (lexically)...
        # and secondarily on the k-mer itself
        self.kmers.sort(key=itemgetter(1, 0))
        with open("minimizers.out", "w") as outfile:
            for match in self.matches:
                outfile.write(">" + str(match.id) + "\n")
                for kmer in match.kmers:
                    outfile.write(str(kmer[0]) + " " + str(kmer[1]) + " " + str(int(kmer[2])))

    # kmer libraries are stored with the name of each repeat prepended with '>',
    # with each kmer-offset pair below on its own line:
    # >repeat1
    # atggtgtggtgatc 0
    # cggcgcatgcgcgc 3
    # >repeat2
    # attgttgttctagc 1
    # aggcgttgcacgat 2
    def parseKmerLibrary(self, filename):
        with open(filename, "r") as kmerFile:
            for line in kmerFile:
                if line[0] == '>':
                    name = line[1:]
                else:
                    kmer, minOffset, isRevComp = line.split()
                    self.kmers[name].append((kmer, int(minOffset), isRevComp))

    ### TEMPORARILY DEPRECATED FOR SPEED
    # returns a k-mer generator
    #def kmerGen(self):
        # only first ten used!
    #    for repeat in self.repeats:
    #        seq = self.chromosomes[repeat.refid][repeat.chromStart-1:repeat.chromEnd]
    #        for i in xrange(len(seq) - self.k + 1):
    #            yield seq[i:i + self.k]

    @staticmethod
    def reverseComplement(seq):
        key = {
            ord('a'): ord('t'), ord('A'): ord('T'), ord('t'): ord('a'), ord('T'): ord('A'),
            ord('c'): ord('g'), ord('C'): ord('G'), ord('g'): ord('c'), ord('G'): ord('C')
            }
        return bytearray([key[seq[len(seq)-i-1]] for i in range(len(seq))])

    # store the indexes of the m-long minimizers of each, using generator list comprehension for efficiency
    # could be done with a nifty one-line generator, but hand-coded for speed
    def minimize(self):
        for match in self.matches:
            # !!! this is dubious - it's unclear whether the sole match's sole sequence identifier represents a file (chrom) or a sequence (seq)
            # for our test case, dm3, the two are the same
            seq = self.chroms[match.chromName].seqs[match.chromName]
            start = match.chromStart - 1
            end = match.chromEnd
            # loop through each of the match's kmers
            for kStart in xrange(start, end - self.k + 1):
                kEnd = kStart + self.k
                minOffset = 0
                isRevComp = False
                currMin = seq[kStart:kStart+self.k]
                # loop through each of the kmer's minimizers
                for mStart in xrange(kStart, kEnd - self.m + 1):
                    testMin = seq[mStart:mStart+self.m]
                    if testMin < currMin:
                        minOffset = mStart - kStart
                        currMin = testMin
                        isRevComp = False
                    revCompMin = self.reverseComplement(seq[mStart:mStart+self.m])
                    if revCompMin < currMin:
                        minOffset = mStart - kStart
                        currMin = revCompMin
                        isRevComp = True
                # self.kmers format: kmer offset, minimizer offset, minimizer is in reverse complement
                match.kmers.append((kStart, minOffset, isRevComp))


def printTree(node, indent=0):
    print ("\t" * indent) + node.cat[-1]
    for child in node.children:
        printTree(child, indent+1)


if __name__ == "__main__":

    chromFilenames = map(lambda x: "dm3/" + x, ["chr2L.fa", "chr2LHet.fa", "chr2R.fa", "chr2RHet.fa", "chr3L.fa", "chr3LHet.fa", "chr3R.fa", "chr3RHet.fa", "chr4.fa", "chrM.fa", "chrUextra.fa", "chrU.fa", "chrX.fa", "chrXHet.fa", "chrYHet.fa"])
    genome = Genome(31, 13, chromFilenames, "dm3/dm3.fa.out")
    printTree(genome.tree)
