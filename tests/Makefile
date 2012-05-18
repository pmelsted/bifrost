CC = g++
CXX = g++
INCLUDES = -I.
CXXFLAGS = -std=c++0x -Wall -Wextra -pedantic $(INCLUDES)


OBJECTS = ../Kmer.o ../hash.o

all: CXXFLAGS += -O3
all: target

target: KmerTest KmerTestExtended

debug: CXXFLAGS += -g -O0
debug: LDFLAGS += -g
debug: target

KmerTest: KmerTest.o 
	$(CC) $(INCLUDES) $(OBJECTS) KmerTest.o -o KmerTest

KmerTestExtended: KmerTestExtended.o 
	$(CC) $(INCLUDES) $(OBJECTS) KmerTestExtended.o -o KmerTestExtended

KmerTest.o: ../Kmer.o ../hash.o
KmerTestExtended.o: ../Kmer.o ../hash.o
Kmer.o: ../Kmer.hpp ../Kmer.cpp
hash.o: ../hash.hpp ../hash.cpp


clean:
	rm -rf *.o
	rm -rf KmerTest KmerTestExtended
