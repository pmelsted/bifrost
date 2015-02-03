# This affects the memory usage of the program
# we use 1 byte for every 4 bp in kmers. Ideally
# this parameter should be a multiple of 4.
# Actual maximum kmer size is 1 less.
MAX_KMER_SIZE = 32

#CC = g++
#CXX = g++
WARNINGS = -Wall -Wno-reorder -Wno-unused-function
INCLUDES = -I. 
CXXFLAGS = -c -std=c++11  $(INCLUDES) -DMAX_KMER_SIZE=$(MAX_KMER_SIZE) $(WARNINGS)
LDFLAGS = -lstdc++ 
LDLIBS  = -lm -lz #-lgomp

all: CXXFLAGS += -O3 
all: target


debug: CXXFLAGS += -g -ggdb  -O0
debug: LDFLAGS += -g -ggdb
debug: target


profile: CXXFLAGS += -p -g -O2
profile: LDFLAGS += -p -g
profile: clean
profile: target

target: BFGraph


OBJECTS = Kmer.o KmerIterator.o  hash.o fastq.o FilterReads.o BuildContigs.o  \
		  CompressedSequence.o Contig.o CompressedCoverage.o  ContigMapper.o

debugtest: CXXFLAGS += -g -O0
debugtest: LDFLAGS += -g
debugtest: debugtest.o $(OBJECTS)
	$(CC) $(INCLUDES) $(OBJECTS) debugtest.o $(LDFLAGS) $(LDLIBS) -o debugtest

BFGraph: BFGraph.o $(OBJECTS)
	$(CXX) $(INCLUDES) $(OBJECTS) BFGraph.o $(LDFLAGS) $(LDLIBS) -o BFGraph -lstdc++



BFGraph.o: BFGraph.cpp
FilterReads.o: FilterReads.cpp  BlockedBloomFilter.hpp
BuildContigs.o: BuildContigs.cpp  BlockedBloomFilter.hpp
fastq.o: fastq.hpp fastq.cpp 
kmer.o: kmer.hpp kmer.cpp
KmerIterator.o: KmerIterator.hpp KmerIterator.cpp
CompressedSequence.o: CompressedSequence.cpp CompressedSequence.hpp
CompressedCoverage.o: CompressedCoverage.cpp CompressedCoverage.hpp
Contig.o: Contig.cpp Contig.hpp 
ContigMapper.o: ContigMapper.cpp ContigMapper.hpp
debugtest.o: debugtest.cpp
hash.o: hash.hpp hash.cpp	


clean:
	rm -f *.o BFGraph
