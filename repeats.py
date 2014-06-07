# WARNING: repeats[:10] used in Genome.__init__() for debugging

from pprint import pprint
from collections import defaultdict


class Repeat:

    def __init__(self, line):
        # parse fields
        (self.swsc, self.pctdiv, self.pctdel, self.pctins, self.refid,
         self.ref_i, self.ref_f, self.ref_remain, self.orient, self.rep_nm,
         self.rep_cl, self.rep_prior, self.rep_i, self.rep_f, self.unk) = line.split()
        # int-ize the reference coordinates
        self.ref_i, self.ref_f = int(self.ref_i), int(self.ref_f)
        # split on dash to give category heirarchy
        self.rep_cl = self.rep_cl.split('/')
        
        # initialize minimizers list, to make it clear we'll be populating it soon
        self.minimizers = None

    def parseSeq(self, seq):
        pass


def parse_fasta(fns):
    ret = defaultdict(bytearray)
    for fn in fns:
        with open(fn) as fh:
            for line in fh:
                if line[0] == '>':
                    name = line[1:].rstrip()
                else:
                    ret[name].extend(line.rstrip())
    return ret


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
        self.categoryNodes = None
        self.chromosomes = defaultdict(bytearray)
        # length of each k-mer
        self.k = k
        # length of each minimizer
        self.m = m

        for filename in chromFilenames:
            with open(filename, "r") as chromFile:
                for line in chromFile:
                    if line[0] == '>':
                        chromName = line[1:].rstrip()
                    else:
                        self.chromosomes[chromName] += line.rstrip()

        for repeat in self.repeats[:10]:
            kmerGen = self.makeKmerGen(self.chromosomes[repeat.refid][repeat.ref_i-1:repeat.ref_f])
            repeat.minimizers = [self.minimize(kmer) for kmer in kmerGen]
            print repeat.minimizers

    # returns a k-mer generator
    def makeKmerGen(self, seq):
        for i in xrange(len(seq) - self.k + 1):
            yield seq[i:i + self.k]

    # returns m-long minimizer of k-mer, using generator list comprehension for efficiency
    def minimize(self, kmer):
        return min((kmer[i:i + self.m] for i in xrange(self.k - self.m + 1)))

    # node definition for below tree
    class CategoryNode:
        def __init__(self, cat):
            self.name = cat[-1]
            self.parentName = "root" if len(cat) < 2 else cat[-2]
            # makeTree eventually puts a pointer here
            self.parent = None
            # makeTree eventually populates this
            self.children = set()

    # each node is a tuple containing the repeat category's name and its associated k-mers
    def makeTree(self):
        CategoryNode = self.CategoryNode    # for brevity
        self.tree = CategoryNode([None, "root"])
        self.categoryNodes = {"root": self.tree}
        # for every unique repeat category
        categories = set((lambda r: r.ref_cl for repeat in self.repeats))
        print type(self.repeats[0].ref_cl)
        for cat in categories:
            CategoryNode = CategoryNode(cat)
            self.categoryNodes[cat.name] = cat
            cat.parent.children.add(cat)


def printTree(node, indent=0):
    print ("\t" * indent), node.name
    for child in root.children:
        printTree(indent + 1, child)


if __name__ == "__main__":

    #types = set(["/".join(repeat.rep_cl) for repeat in repeats])
    #pprint(types)

    chromFilenames = ["chr2L.fa", "chr2LHet.fa", "chr2R.fa", "chr2RHet.fa", "chr3L.fa", "chr3LHet.fa", "chr3R.fa", "chr3RHet.fa", "chr4.fa", "chrM.fa", "chromFa.tar", "chrUextra.fa", "chrU.fa", "chrX.fa", "chrXHet.fa", "chrYHet.fa"]
    genome = Genome(chromFilenames, "dm3.fa.out", 31, 13)
    genome.makeTree()
    printTree(genome.tree)

    repeatFile.close()
