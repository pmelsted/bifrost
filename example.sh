#!/bin/bash
set -e  # Exit on errors

#0) Make sure the program has been compiled
if [! -f "BFGraph" ]
then
    echo "BFGraph not found, did you run make?"
    exit 1
fi


#1) Filter the reads in 1 thread with kmersize 31
./BFGraph filter example/tinyread_*.fq -k 31 -t 1 -o example/output/tiny.bf -n 8000 -N 4000 -v

#2) Make the contigs, kmersize must be the same as in #1
./BFGraph contigs example/tinyread_*.fq -k 31 -t 1 -f example/output/tiny.bf -o example/output/tiny -v

#3) Create the graph
./make_graph.py example/output/tiny

#4) Tell the user how to create a PNG from the .dot file
echo -e "\nYou can now run: \`./dot2png.sh example/output/tiny\`\nto create 'example/output/tiny.png', a PNG of the De Brujin graph"
