import sys
import matplotlib
from collections import defaultdict

if len(sys.argv) != 2:
    print "must supply genome to be analyzed"
    sys.exit(1)

genomeName = sys.argv[1]

lcaFreq = defaultdict(int)
with open(genomeName + ".mins", "r") as f:
    for line in f:
        if line[:2] == "\t":
            fields = line.split()
            if len(fields) != 2:
                print "ERROR: malformed line encountered:", line
                sys.exit(1)
            kmer, lca = fields[0], fields[1]
            lcaFreq[lca] += 1

with open(genomeName + ".lcafreq", "w") as out:
    out.writelines([lca + ' ' + str(count) + '\n' for lca, count in lcaFreq.iteritems()])


matchSizes = defaultdict(int)
nodeSizes = defaultdict(int)
with open(genomeName + '/' + genomeName + ".fa.out", "r") as f:
    _ = [f.readline() for _ in range(3)]
    for line in f:
        fields = line.split()
        if len(fields) != 15:
            print "malformed match line"
            sys.exit(1)
        seqStart, seqEnd = int(fields[5]) - 1, int(fields[6])
        matchSizes[seqEnd - seqStart] += 1

with open(genomeName + ".matchSizes", "w") as f:
    out.writelines((str(size) + ' ' + str(count) + '\n' for size, count in matchSizes.iteritems()))


