#!/bin/bash

echo "Running KmerTest"
./KmerTest 10
echo -e "\n"

echo "Running KmerTestExtended"
./KmerTestExtended 31 5
echo -e "\n"

echo "Running CompressedSequenceTest"
./CompressedSequenceTest 20
echo -e "\n"

echo "Running BloomFilterTest"
./BloomFilterTest 4
echo -e "\n"

echo "Running BlockedBloomFilterTest"
./BlockedBloomFilterTest 4
echo -e "\n"

echo "Running KmerMapperTest"
./KmerMapperTest
