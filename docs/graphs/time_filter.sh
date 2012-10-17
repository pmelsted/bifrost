#!/bin/bash
set -e

# Must be run inside directory docs/graphs
cd ../../

if [ ! -f "data/small.fa" ] 
then
    echo "missing file ../../data/small.fa"
    cd docs/graphs
    exit 1
fi


for threads in {1..8}; do \
    beginning=`date +"%s.%N"`
    echo -ne "$threads "
    ./BFGraph filter -k 31 -t $threads -o data/big.bf -n 8000000 -N 4000000 data/small.fa 2>/dev/null
    end=`date +"%s.%N"`
    echo `python -c "print '%.3f' % ($end-$beginning),"`
done

cd docs/graphs
