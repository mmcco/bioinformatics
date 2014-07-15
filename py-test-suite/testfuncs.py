import itertools

def revComp(seq):
    seq = bytearray(seq)
    key = {
        ord('a'): ord('t'),
        ord('t'): ord('a'),
        ord('c'): ord('g'),
        ord('g'): ord('c'),
        ord('n'): ord('n'),
        }
    return bytearray([key[seq[len(seq)-i-1]] for i in range(len(seq))])

# store the indexes of the m-long minimizers of each, using generator list comprehension for efficiency
# could be done with a nifty one-line generator, but hand-coded for speed
def mini(kmer, m):
    if len(kmer) < m:
        raise Exception("error in testfuncs.mini(): len(kmer) < m")
    kmer = bytearray(kmer)
    return min((min(kmer[i:i+m], revComp(kmer[i:i+m])) for i in xrange(len(kmer) - m + 1)))

def miniAll(kmer, m):
    kmer = bytearray(kmer)
    minGen = (min((kmer[i:i+m], i, False), (revComp(kmer[i:i+m]), i, True)) for i in xrange(len(kmer) - m + 1))
    return min(minGen)
