cmake_minimum_required(VERSION 2.8.12)

project(Bifrost)

# Can be replaced with a larger multiple of 32. Actual maximum k-mer size is MAX_KMER_SIZE-1.
#set(MAX_KMER_SIZE "32")
SET(MAX_KMER_SIZE "32" CACHE STRING "MAX_KMER_SIZE")

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c11 -march=native -DXXH_NAMESPACE=BIFROST_HASH_")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -march=native -DXXH_NAMESPACE=BIFROST_HASH_")

#set_property(SOURCE xxhash.c APPEND_STRING PROPERTY COMPILE_FLAGS " -Wall -Wextra -Wcast-qual -Wcast-align -Wshadow -Wstrict-aliasing=1 -Wswitch-enum -Wdeclaration-after-statement -Wstrict-prototypes -Wundef")
set_property(SOURCE BlockedBloomFilter.cpp APPEND_STRING PROPERTY COMPILE_FLAGS " -funroll-loops")

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release")
endif(NOT CMAKE_BUILD_TYPE)

if(CMAKE_BUILD_TYPE MATCHES Debug)
	message("Debug mode")
	add_compile_options(-g)
else(CMAKE_BUILD_TYPE MATCHES Debug)
	if(CMAKE_BUILD_TYPE MATCHES Profile)
		message("Profile mode")
		add_compile_options(-pg)
		set(CMAKE_SHARED_LINKER_FLAGS "-pg")
		set(CMAKE_EXE_LINKER_FLAGS "-pg")
	else(CMAKE_BUILD_TYPE MATCHES Profile)
		message("Release mode")
		add_compile_options(-O3)
	endif(CMAKE_BUILD_TYPE MATCHES Profile)
endif(CMAKE_BUILD_TYPE MATCHES Debug)

add_subdirectory(src)
