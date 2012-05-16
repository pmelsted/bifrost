CC = g++
INCLUDES = -I.


OBJECTS = ../Kmer.o ../hash.o

all: CXXFLAGS += -O3
all: target

target: KmerTest

debug: CXXFLAGS += -g -O0
debug: LDFLAGS += -g
debug: target

KmerTest: KmerTest.o 
	$(CC) $(INCLUDES) $(OBJECTS) KmerTest.o -o KmerTest

KmerTest.o: ../Kmer.o ../hash.o
Kmer.o: ../Kmer.hpp ../Kmer.cpp
hash.o: ../hash.hpp ../hash.cpp


clean:
	rm -rf *.o
	rm -rf KmerTest
