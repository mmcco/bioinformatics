'''
A basic Bloom filter implementation in order to
teach the concept. It will soon be incorporated
into an iPython notebook.
'''

class BloomFilter:

    def __init__(self, m, k):
        self.m = m
        self.k = k
        self.bitArray = [0 for _ in range(0, m)]

    def getIndices(self, val):
        indices = []
        for i in range(0, self.k):
            thisHash = hash(str(i) + val)
            index = thisHash % self.m
            indices.append(index)
        return indices

    def add(self, val):
        indices = self.getIndices(val)
        for index in indices:
            self.bitArray[index] = 1

    def check(self, val):
        indices = self.getIndices(val)
        for index in indices:
            if self.bitArray[index]  == 0:
                return False
        return True


test = BloomFilter(256, 3)
print test.bitArray
print test.getIndices("testing")
test.add("testing")
print test.bitArray
print test.check("testing")
print test.check("not in the set")
