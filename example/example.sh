#!/bin/bash
set -e  # Exit on errors

#0) Check that make has been run
if [ ! -f "../BFGraph" ]
then
    echo "../BFGraph not found, did you run make in the parent directory?"
    exit 1
fi

if [ ! -d "output" ]
then
    echo "Creating directory: output"
    mkdir output
fi


#1) Filter the reads in 1 thread with kmersize 31
../BFGraph filter -k 31 -t 1 -N 8000000 -n 4000000 -o output/tiny.bf tinyread_*.fq -v

#2) Make the contigs, kmersize must be the same as in #1
../BFGraph contigs -t 1 -k 31 -f output/tiny.bf -o output/tiny tinyread_*.fq -v

#3) Create the graph
../make_graph.py output/tiny

#4) Get the contig strings from the contig file without >contigXX lines
grep -v contig output/tiny.contigs > output/tiny.contigs.cleaned

#5) Sort the reps of the contigs alphabetically into the file data/tiny.contigs.sorted
python ../sort_contigs.py output/tiny.contigs.cleaned output/tiny.contigs.sorted

#6) Diff the sorted contigs against the correct output
echo "There should be no diff. Diff start:"
diff output/tiny.contigs.sorted tinycontigs.correct
echo "End of diff"

#7) Tell the user how to create a PNG from the .dot file
echo "You can now run: '../dot2png.sh output/tiny' to create 'output/tiny.png', a PNG of the De Brujin graph"
