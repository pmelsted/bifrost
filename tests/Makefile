MAX_KMER_SIZE=96

CC = g++
CXX = g++
INCLUDES = -I.
CXXFLAGS = -c -Wall -Wno-reorder $(INCLUDES) -fPIC -DMAX_KMER_SIZE=$(MAX_KMER_SIZE)
LDFLAGS = -lz

EXECUTABLES = KmerTest KmerTestExtended CompressedSequenceTest BloomFilterTest KmerMapperTest KmerIteratorTest \
			  CompressedCoverageTest KmerHashTest

all: CXXFLAGS += -O3
all: target

target: $(EXECUTABLES) 

debug: CXXFLAGS += -g -O0
debug: LDFLAGS += -g
debug: target

profile: CXXFLAGS += -p -g -O2
profile: LDFLAGS += -p -g
profile: clean
profile: target

OBJECTS = ../FindContig.o ../KmerMapper.o ../Kmer.o ../hash.o ../CompressedSequence.o ../Contig.o  ../KmerIterator.o ../CompressedCoverage.o ../fastq.o

KmerTest: KmerTest.o $(OBJECTS)
	$(CC) $(INCLUDES) $(OBJECTS) KmerTest.o $(LDFLAGS) -o KmerTest

KmerHashTest: KmerHashTest.o $(OBJECTS)
	$(CC) $(INCLUDES) $(OBJECTS) KmerHashTest.o $(LDFLAGS) -o KmerHashTest -lstdc++


KmerTestExtended: KmerTestExtended.o $(OBJECTS)
	$(CC) $(INCLUDES) $(OBJECTS) KmerTestExtended.o $(LDFLAGS) -o KmerTestExtended

CompressedSequenceTest: CompressedSequenceTest.o $(OBJECTS)
	$(CC) $(INCLUDES) $(OBJECTS) CompressedSequenceTest.o $(LDFLAGS) -o CompressedSequenceTest

BloomFilterTest: BloomFilterTest.o $(OBJECTS)
	$(CC) $(INCLUDES) $(OBJECTS) BloomFilterTest.o $(LDFLAGS) -o BloomFilterTest

BlockedBloomFilterTest: BlockedBloomFilterTest.o $(OBJECTS)
	$(CC) $(INCLUDES) $(OBJECTS) BlockedBloomFilterTest.o $(LDFLAGS) -o BlockedBloomFilterTest

KmerMapperTest: KmerMapperTest.o $(OBJECTS)
	$(CC) $(INCLUDES) $(OBJECTS) KmerMapperTest.o $(LDFLAGS) -o KmerMapperTest

KmerIteratorTest: KmerIteratorTest.o $(OBJECTS)
	$(CC) $(INCLUDES) $(OBJECTS) KmerIteratorTest.o $(LDFLAGS) -o KmerIteratorTest

CompressedCoverageTest: CompressedCoverageTest.o $(OBJECTS)
	$(CC) $(INCLUDES) $(OBJECTS) CompressedCoverageTest.o $(LDFLAGS) -o CompressedCoverageTest


KmerTest.o: ../Kmer.o ../hash.o
KmerHashTest.o: ../Kmer.o ../hash.o ../KmerIterator.o
KmerTestExtended.o: ../Kmer.o ../hash.o
CompressedSequenceTest.o: ../CompressedSequence.o
BloomFilterTest.o: ../BloomFilter.o 
BlockedBloomFilterTest.o: ../BlockedBloomFilter.o
KmerMapperTest.o: ../KmerMapper.o
KmerIteratorTest.o: ../KmerIterator.o
CompressedCoverageTest.o: ../CompressedCoverage.o


../Kmer.o: ../Kmer.hpp ../Kmer.cpp
../hash.o: ../hash.hpp ../hash.cpp
../CompressedSequence.o: ../Common.hpp ../CompressedSequence.hpp ../CompressedSequence.cpp
../BloomFilter.o: ../BloomFilter.hpp
../BlockedBloomFilter.o: ../BlockedBloomFilter.hpp
../BlockedBloomFilter.hpp: ../BloomFilter.hpp
../KmerMapper.o: ../KmerMapper.hpp ../KmerMapper.cpp
../KmerIterator.o: ../KmerIterator.hpp ../KmerIterator.cpp
../CompressedCoverage.o: ../CompressedCoverage.hpp ../CompressedCoverage.cpp
../Contig.o: ../Contig.hpp ../Contig.cpp
../FindContig.o: ../FindContig.hpp ../FindContig.cpp

clean:
	rm -f *.o ../*.o $(EXECUTABLES) 
