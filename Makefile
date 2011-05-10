CC = g++
CXX = g++
INCLUDES = -I.
CXXFLAGS = -c -Wall -Wno-reorder $(INCLUDES)
LDFLAGS =
LDLIBS  = -lm -lz




all: CXXFLAGS += -O3
all: target


debug: CXXFLAGS += -g -O0
debug: LDFLAGS += -g
debug: target

profile: CXXFLAGS += -pg -O2
profile: LDFLAGS += -pg
profile: clean
profile: target

naive: Naive

target: BFCounter

OBJECTS =  kmer.o hash.o bloom_filter.o fastq.o 

testread: testread.o $(OBJECTS)
	$(CC) $(INCLUDES) $(LDFLAGS) $(LDLIBS) $(OBJECTS) testread.o -o testread

BFCounter: BFCounter.o $(OBJECTS)
	$(CC) $(INCLUDES) $(LDFLAGS) $(LDLIBS) $(OBJECTS) BFCounter.o -o BFCounter

Naive: Naive.o $(OBJECTS)
	$(CC) $(INCLUDES) $(LDFLAGS) $(LDLIBS) $(OBJECTS) Naive.o -o Naive


BFCounter.o: BFCounter.cpp
fastq.o: fastq.hpp fastq.cpp 
kmer.o: kmer.hpp kmer.cpp
bloom_filter.o: bloom_filter.hpp bloom_filter.cpp
hash.o: hash.hpp hash.cpp	
Naive.o: Naive.cpp

clean:
	rm -rf *.o
	rm -rf BFCounter
	rm -rf Naive
