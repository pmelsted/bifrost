cmake_minimum_required(VERSION 2.8.12)

project(Bifrost)

# To enable a larger default k-mer size, replace MAX_KMER_SIZE with a larger multiple of 32: actual maximum k-mer size will be MAX_KMER_SIZE-1.
SET(MAX_KMER_SIZE "32" CACHE STRING "MAX_KMER_SIZE")
# Enable architecture optimizations
SET(ENABLE_NATIVE_ARCH "ON" CACHE STRING "ENABLE_NATIVE_ARCH")
# Enable AVX2 instructions
SET(ENABLE_AVX2 "ON" CACHE STRING "ENABLE_AVX2")

# Set some default compile flags
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c11 -DXXH_NAMESPACE=BIFROST_HASH_")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -DXXH_NAMESPACE=BIFROST_HASH_")

set_property(SOURCE BlockedBloomFilter.cpp APPEND_STRING PROPERTY COMPILE_FLAGS " -funroll-loops")

if(ENABLE_NATIVE_ARCH MATCHES "ON")
	set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -march=native")
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native")
else(ENABLE_NATIVE_ARCH MATCHES "ON")
	message("Disabling native architecture compilation (including AVX2)")
endif(ENABLE_NATIVE_ARCH MATCHES "ON")

if(ENABLE_AVX2 MATCHES "OFF")
	message("Disabling AVX2 instructions")
	set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -mno-avx2")
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mno-avx2")
endif(ENABLE_AVX2 MATCHES "OFF")

# Manages build types 
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release")
endif(NOT CMAKE_BUILD_TYPE)

if(CMAKE_BUILD_TYPE MATCHES Debug)
	message("Build type: Debug")
	add_compile_options(-g)
else(CMAKE_BUILD_TYPE MATCHES Debug)
	if(CMAKE_BUILD_TYPE MATCHES Profile)
		message("Build type: Profiling")
		add_compile_options(-pg)
		set(CMAKE_SHARED_LINKER_FLAGS "-pg")
		set(CMAKE_EXE_LINKER_FLAGS "-pg")
	else(CMAKE_BUILD_TYPE MATCHES Profile)
		message("Build type: Release")
		add_compile_options(-O3)
	endif(CMAKE_BUILD_TYPE MATCHES Profile)
endif(CMAKE_BUILD_TYPE MATCHES Debug)

MATH(EXPR PRINT_MAX_KMER_SIZE "${MAX_KMER_SIZE}-1")
message("Maximum k-mer size is " ${PRINT_MAX_KMER_SIZE})

add_subdirectory(src)
