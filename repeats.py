from pprint import pprint
from collections import defaultdict
from operator import itemgetter


class Repeat:

    def __init__(self, line):
        # parse fields
        (self.swsc, self.pctdiv, self.pctdel, self.pctins, self.refid,
         self.ref_i, self.ref_f, self.ref_remain, self.orient, self.rep_nm,
         self.rep_cl, self.rep_prior, self.rep_i, self.rep_f, self.unk) = line.split()
        # int-ize the reference coordinates
        self.ref_i, self.ref_f = int(self.ref_i), int(self.ref_f)


class Genome:

    def __init__(self, chromFilenames, repeatFilename, k, m):
        with open(repeatFilename, "r") as repeatFile:
            repeatFile.readline()
            repeatFile.readline()
            repeatFile.readline()
            self.repeats = [Repeat(line) for line in repeatFile]
        # contains a pointer to the tree root, populated by makeTree()
        self.tree = None
        # a dict of CategoryNodes indexed by name, populated by makeTree()
        self.catNodes = None
        self.chromosomes = defaultdict(bytearray)
        # length of each k-mer
        self.k = k
        # length of each minimizer
        self.m = m
        # populated later by minimize()
        self.kmers = []

        for filename in chromFilenames:
            with open(filename, "r") as chromFile:
                for line in chromFile:
                    if line[0] == '>':
                        chromName = line[1:].rstrip()
                    else:
                        self.chromosomes[chromName] += line.rstrip()

        # kmers stored in a tuple: (kmer, minimizer offset, chromName)
        self.minimize()
        # use Kraken technique, sorting primarily on minimizer (lexically)...
        # and secondarily on the k-mer itself
        self.kmers.sort(key=itemgetter(1, 0))
        with open("minimizers.out", "w") as outfile:
            outfile.writelines((t[0] + ' ' + str(t[1]) for t in self.kmers))

    # returns a k-mer generator
    def kmerGen(self):
        # only first ten used!
        for repeat in self.repeats[:5000]:
            seq = self.chromosomes[repeat.refid][repeat.ref_i-1:repeat.ref_f]
            for i in xrange(len(seq) - self.k + 1):
                yield seq[i:i + self.k]

    # returns the index of m-long minimizer of k-mer, using generator list comprehension for efficiency
    # could be done with a nifty one-line generator, but hand-coded for speed
    def minimize(self):
        temp = 0    # DELETE ME
        for repeat in self.repeats:
            temp += 1    # DELETE ME
            chrom = self.chromosomes[repeat.refid]
            start = repeat.ref_i - 1
            end = repeat.ref_f
            # loop through each of the repeat's kmers
            for kStart in xrange(start, end - self.k + 1):
                kEnd = kStart + self.k
                min_offset = 0
                curr_min = chrom[start:start+self.k]
                # loop through each of the kmer's minimizers
                for mStart in xrange(kStart, kEnd - self.m + 1):
                    if chrom[mStart:mStart+self.m] < chrom[kStart+min_offset:kStart+min_offset+self.m]:
                        min_offset = mStart - kStart
                self.kmers.append((chrom[kStart:kEnd], min_offset, repeat.refid))
                if temp % 1000 == 0 and kStart == start:
                    print temp

    # node definition for below tree
    class CategoryNode:
        def __init__(self, cat):
            # nameList and name with the name as a list (implyig heirarchy) and as a string respectively
            # a bit of a hack, but constantly joining and splitting takes time, and there aren't many cats
            self.nameList = cat
            self.name = '/'.join(cat)
            self.parentName = None if self.name == "root" else "root" if len(cat) < 2 else cat[-2]
            # makeTree populates self.parent and self.children
            self.parent = None
            self.children = set()

    # each node is a tuple containing the repeat category's name and its associated k-mers
    def makeTree(self):
        CategoryNode = self.CategoryNode    # for brevity
        self.tree = CategoryNode(["root"])
        self.catNodes = {"root": self.tree}
        otherID = 0    # a counter for handling the "Other" category below
        # generator of all repeat categories
        categories = (repeat.rep_cl for repeat in self.repeats)
        for cats in categories:
            # split cats into heirarchy list
            catList = cats.split('/')
            # must add every heirarchy level to the tree, hence the branch loop
            branch = (catList[:i+1] for i in xrange(len(catList)))
            for cat in branch:
                # must remove duplicates manually because Python doesn't allow sets of lists
                if '/'.join(cat) in self.catNodes:
                    continue
                else:
                    # an ugly conditional to ensure that unknowns are not grouped as a category,
                    # but are rather children of root
                    if cat[-1] == "Other":
                        cat[-1] = "Unknown-" + str(otherID)
                        otherID += 1
                    catNode = CategoryNode(cat)
                    self.catNodes[catNode.name] = catNode
        # now loop again to populate parents and children
        for catNode in self.catNodes.values():
            if catNode.parentName:
                catNode.parent = self.catNodes[catNode.parentName]
                catNode.parent.children.add(catNode)
            else:
                catNode.parent = None


def printTree(node, indent=0):
    print ("\t" * indent) + node.nameList[-1]
    for child in node.children:
        printTree(child, indent+1)


if __name__ == "__main__":

    chromFilenames = ["chr2L.fa", "chr2LHet.fa", "chr2R.fa", "chr2RHet.fa", "chr3L.fa", "chr3LHet.fa", "chr3R.fa", "chr3RHet.fa", "chr4.fa", "chrM.fa", "chrUextra.fa", "chrU.fa", "chrX.fa", "chrXHet.fa", "chrYHet.fa"]
    genome = Genome(chromFilenames, "dm3.fa.out", 31, 13)
    genome.makeTree()
    #for kmer in genome.kmers[:100]:
    #    print kmer[0], "    ---->    ", kmer[1]
    #printTree(genome.tree)
