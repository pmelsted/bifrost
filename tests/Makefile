CC = g++
INCLUDES = -I.


OBJECTS = ../Kmer.o ../hash.o

KmerTest: KmerTest.o 
	$(CC) $(INCLUDES) $(OBJECTS) KmerTest.o -o KmerTest

KmerTest.o: ../Kmer.o ../hash.o
Kmer.o: ../Kmer.hpp ../Kmer.cpp
hash.o: ../hash.hpp ../hash.cpp

