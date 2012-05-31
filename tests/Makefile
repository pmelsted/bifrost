CC = g++
CXX = g++
INCLUDES = -I.
CXXFLAGS = -pg -Wall -Wextra -pedantic $(INCLUDES)

EXECUTABLES = KmerTest KmerTestExtended CompressedSequenceTest BloomFilterTest BlockedBloomFilterTest KmerMapperTest

all: CXXFLAGS += -O3
all: target

target: $(EXECUTABLES) 

debug: CXXFLAGS += -g -O0
debug: LDFLAGS += -g
debug: target

OBJECTS = ../Kmer.o ../hash.o ../CompressedSequence.o ../Contig.o ../KmerMapper.o

KmerTest: KmerTest.o $(OBJECTS)
	$(CC) $(INCLUDES) $(OBJECTS) KmerTest.o -o KmerTest

KmerTestExtended: KmerTestExtended.o $(OBJECTS)
	$(CC) $(INCLUDES) $(OBJECTS) KmerTestExtended.o -o KmerTestExtended

CompressedSequenceTest: CompressedSequenceTest.o $(OBJECTS)
	$(CC) $(INCLUDES) $(OBJECTS) CompressedSequenceTest.o -o CompressedSequenceTest

BloomFilterTest: BloomFilterTest.o $(OBJECTS)
	$(CC) $(INCLUDES) $(OBJECTS) BloomFilterTest.o -o BloomFilterTest

BlockedBloomFilterTest: BlockedBloomFilterTest.o $(OBJECTS)
	$(CC) $(INCLUDES) $(OBJECTS) BlockedBloomFilterTest.o -o BlockedBloomFilterTest

KmerMapperTest: KmerMapperTest.o $(OBJECTS)
	$(CC) $(INCLUDES) $(OBJECTS) KmerMapperTest.o -o KmerMapperTest

KmerTest.o: ../Kmer.o ../hash.o
KmerTestExtended.o: ../Kmer.o ../hash.o
CompressedSequenceTest.o: ../CompressedSequence.o
BloomFilterTest.o: ../BloomFilter.o 
BlockedBloomFilterTest.o: ../BlockedBloomFilter.o
KmerMapperTest.o: ../KmerMapper.o

../Kmer.o: ../Kmer.hpp ../Kmer.cpp
../hash.o: ../hash.hpp ../hash.cpp
../CompressedSequence.o: ../Common.hpp ../CompressedSequence.hpp ../CompressedSequence.cpp
../BloomFilter.o: ../BloomFilter.hpp
../BlockedBloomFilter.o: ../BlockedBloomFilter.hpp
../KmerMapper.o: ../KmerMapper.hpp ../KmerMapper.cpp


clean:
	rm -rf *.o
	rm -rf $(EXECUTABLES) 
