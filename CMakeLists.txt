cmake_minimum_required(VERSION 2.8.12)

project(BFGraph)

# This affects the memory usage of the program
# we use 1 byte for every 4 bp in kmers. Ideally
# this parameter should be a multiple of 4.
# Actual maximum kmer size is 1 less.
set(MAX_KMER_SIZE "32")

SET(CMAKE_C_FLAGS "-std=c99 -DXXH_NAMESPACE=BFGRAPH_HASH_")
SET(CMAKE_CXX_FLAGS "-std=c++11 -DXXH_NAMESPACE=BFGRAPH_HASH_")

set_property(SOURCE xxhash.c APPEND_STRING PROPERTY COMPILE_FLAGS " -Wall -Wextra -Wcast-qual -Wcast-align -Wshadow -Wstrict-aliasing=1 -Wswitch-enum -Wdeclaration-after-statement -Wstrict-prototypes -Wundef")

set_property(SOURCE BlockedBloomFilter.hpp APPEND_STRING PROPERTY COMPILE_FLAGS " -funroll-loops")

if(CMAKE_BUILD_TYPE MATCHES Debug)
	message("debug mode")
	add_compile_options(-g)
else(CMAKE_BUILD_TYPE MATCHES Debug)
	if(CMAKE_BUILD_TYPE MATCHES Profile)
		message("profile mode")
		add_compile_options(-pg)
		set (CMAKE_SHARED_LINKER_FLAGS "-pg")
		set (CMAKE_EXE_LINKER_FLAGS "-pg")
	else(CMAKE_BUILD_TYPE MATCHES Profile)
		message("release mode")
		add_compile_options(-O2)
	endif(CMAKE_BUILD_TYPE MATCHES Profile)
endif(CMAKE_BUILD_TYPE MATCHES Debug)

add_subdirectory(src)
