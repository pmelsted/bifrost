import sys
sys.path.append("../")

import graph  # swig generated

k = 31
s = "ACGTACGTACGTACGTACGTACGTACGTACG"
graph.Kmer.set_k(k)
km = graph.Kmer(s[:k])

bf = graph.BloomFilter()
bf.open('testBloom.bf')
bf.count()

print km
print km in bf
