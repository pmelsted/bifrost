#!/bin/bash
set -e  # Exit on errors

#0) Make sure the program has been compiled
if [ ! -f "BFGraph" ]
then
    echo "BFGraph not found, did you run make?"
    exit 1
fi

#1) Filter the reads with kmersize 31y (this runs in one thread since parameter -t is not given)
./BFGraph filter example/tinyread_*.fq -k 31 -o example/output/tiny.bf -n 8000 -N 4000 -v

#2) Make the contigs, kmersize must be the same as in #1 (this runs in one thread since parameter -t is not given)
./BFGraph contigs example/tinyread_*.fq -k 31 -f example/output/tiny.bf -o example/output/tiny -v
