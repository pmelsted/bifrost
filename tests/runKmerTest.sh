#!/bin/bash

# Use with 1 number as an argument, f.x. ./runKmerTest.sh 10
g++ -o KmerTest KmerTest.cpp ../Kmer.cpp ../hash.cpp
./KmerTest $@
