cmake_minimum_required(VERSION 2.8.12)

project(BFGraph)

# This affects the memory usage of the program
# we use 1 byte for every 4 bp in kmers. Ideally
# this parameter should be a multiple of 4.
# Actual maximum kmer size is 1 less.
set( MAX_KMER_SIZE "32")

add_compile_options(-std=c++11)
IF(CMAKE_BUILD_TYPE MATCHES Debug)
    message("debug mode")
		add_compile_options(-O0 -g)
ELSE(CMAKE_BUILD_TYPE MATCHES Release)
		add_compile_options(-O3)
ENDIF(CMAKE_BUILD_TYPE MATCHES Debug)

if(CMAKE_BUILD_TYPE MATCHES Profile)
	add_compile_options(-g)
endif(CMAKE_BUILD_TYPE MATCHES Profile)

add_subdirectory(src)
