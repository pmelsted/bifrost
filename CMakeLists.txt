cmake_minimum_required(VERSION 3.13)
set(CMAKE_CONFIGURATION_TYPES "${CMAKE_CONFIGURATION_TYPES};Profile")

project(Bifrost)

if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    set(CMAKE_C_FLAGS_PROFILE "${CMAKE_C_FLAGS_RELEASE} -pg")
    set(CMAKE_CXX_FLAGS_PROFILE "${CMAKE_CXX_FLAGS_RELEASE} -pg")
    set(CMAKE_SHARED_LINKER_FLAGS_PROFILE "${CMAKE_SHARED_LINKER_FLAGS_RELEASE} -pg")
    set(CMAKE_EXE_LINKER_FLAGS_PROFILE "${CMAKE_EXE_LINKER_FLAGS_RELEASE} -pg")
endif()

# Can be replaced with a larger and appropriate number (a power of 2). Actual k-mer size is MAX_KMER_SIZE-1.
set(MAX_KMER_SIZE    "32")

add_subdirectory(src)
