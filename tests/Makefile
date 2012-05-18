CC = g++
CXX = g++
INCLUDES = -I.
CXXFLAGS = -Wall -Wextra -pedantic $(INCLUDES)



all: CXXFLAGS += -O3
all: target

target: KmerTest KmerTestExtended CompressedSequenceTest

debug: CXXFLAGS += -g -O0
debug: LDFLAGS += -g
debug: target

OBJECTS = ../Kmer.o ../hash.o ../CompressedSequence.o

KmerTest: KmerTest.o $(OBJECTS)
	$(CC) $(INCLUDES) $(OBJECTS) KmerTest.o -o KmerTest

KmerTestExtended: KmerTestExtended.o $(OBJECTS)
	$(CC) $(INCLUDES) $(OBJECTS) KmerTestExtended.o -o KmerTestExtended

CompressedSequenceTest: CompressedSequenceTest.o $(OBJECTS)
	$(CC) $(INCLUDES) $(OBJECTS) CompressedSequenceTest.o -o CompressedSequenceTest

KmerTest.o: ../Kmer.o ../hash.o
KmerTestExtended.o: ../Kmer.o ../hash.o
CompressedSequenceTest.o: ../CompressedSequence.o

../Kmer.o: ../Kmer.hpp ../Kmer.cpp
../hash.o: ../hash.hpp ../hash.cpp
../CompressedSequence.o: ../Common.hpp ../CompressedSequence.hpp ../CompressedSequence.cpp


clean:
	rm -rf *.o
	rm -rf KmerTest KmerTestExtended CompressedSequenceTest
