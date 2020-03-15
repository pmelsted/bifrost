file(GLOB sources *.cpp roaring.c)
file(GLOB headers *.h *.hpp *.hh *.tcc)

list(REMOVE_ITEM sources Bifrost.cpp)

add_definitions(-DMAX_KMER_SIZE=${MAX_KMER_SIZE})

if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    # Make sure to add trailing space
    set(PROFILE_FLAGS "$<$<CONFIG:Profile>:-pg> ")
else()
    set(PROFILE_FLAGS "")
endif()

set(BIFROST_COMPILER_FLAGS "${PROFILE_FLAGS}-march=native")
set(BIFROST_LINKER_FLAGS "${PROFILE_FLAGS}")
set(BIFROST_COMPILER_DEFINITIONS "XXH_NAMESPACE=BIFROST_HASH_")

add_library(bifrost_static STATIC ${sources} ${headers})
add_library(bifrost_dynamic SHARED ${sources} ${headers})

set_target_properties(bifrost_static PROPERTIES OUTPUT_NAME "bifrost")
set_target_properties(bifrost_dynamic PROPERTIES OUTPUT_NAME "bifrost")

target_include_directories(bifrost_static PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_include_directories(bifrost_dynamic PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

target_compile_features(bifrost_static PUBLIC cxx_std_11)
target_compile_features(bifrost_dynamic PUBLIC cxx_std_11)

target_compile_options(bifrost_static PUBLIC "${BIFROST_COMPILER_FLAGS}")
target_compile_options(bifrost_dynamic PUBLIC "${BIFROST_COMPILER_FLAGS}")

target_link_options(bifrost_static PUBLIC "${BIFROST_LINKER_FLAGS}")
target_link_options(bifrost_dynamic PUBLIC "${BIFROST_LINKER_FLAGS}")

target_compile_definitions(bifrost_static PRIVATE "${BIFROST_COMPILER_DEFINITIONS}")
target_compile_definitions(bifrost_dynamic PRIVATE "${BIFROST_COMPILER_DEFINITIONS}")

add_executable(Bifrost Bifrost.cpp)

#find_library(JEMALLOC_STATIC libjemalloc_pic.a)
#if (JEMALLOC_STATIC)
#	message("-- Found Jemalloc: TRUE (static)")
#	target_link_libraries(bifrost_static ${JEMALLOC_STATIC} -ldl -lrt)
#	target_link_libraries(bifrost_dynamic ${JEMALLOC_STATIC} -ldl -lrt)
#else()
#	find_library(JEMALLOC_SHARED jemalloc)
#	if (JEMALLOC_SHARED)
#		message("-- Found Jemalloc: TRUE (shared)")
#		target_link_libraries(bifrost_static jemalloc -lrt)
#		target_link_libraries(bifrost_dynamic jemalloc -lrt)
#	else()
#		message("-- Found Jemalloc: FALSE")
#	endif(JEMALLOC_SHARED)
#endif(JEMALLOC_STATIC)

#find_library(CROARING libroaring.a)
#if (CROARING)
#	target_link_libraries(bifrost_static ${CROARING})
#	target_link_libraries(bifrost_dynamic ${CROARING})
#else()
#	target_link_libraries(bifrost_static roaring)
#	target_link_libraries(bifrost_dynamic roaring)
#endif(CROARING)

find_package(Threads REQUIRED)
target_link_libraries(bifrost_static ${CMAKE_THREAD_LIBS_INIT})
target_link_libraries(bifrost_dynamic ${CMAKE_THREAD_LIBS_INIT})

find_package(ZLIB REQUIRED)
target_link_libraries(bifrost_static ${ZLIB_LIBRARIES})
target_link_libraries(bifrost_dynamic ${ZLIB_LIBRARIES})

if (ZLIB_FOUND)
    include_directories(${ZLIB_INCLUDE_DIRS})
else()
    message(FATAL_ERROR "zlib not found. Required for to output files")
endif(ZLIB_FOUND)

target_link_libraries(Bifrost bifrost_dynamic)
target_compile_features(Bifrost PUBLIC cxx_std_11)
target_compile_options(Bifrost PUBLIC "${BIFROST_COMPILER_FLAGS}")
target_link_options(Bifrost PUBLIC "${BIFROST_LINKER_FLAGS}")
target_compile_definitions(Bifrost PRIVATE "${BIFROST_COMPILER_DEFINITIONS}")

set_property(SOURCE BlockedBloomFilter.cpp APPEND_STRING PROPERTY COMPILE_FLAGS " -funroll-loops")

install(TARGETS Bifrost DESTINATION bin)
install(TARGETS bifrost_dynamic DESTINATION lib)
install(TARGETS bifrost_static DESTINATION lib)
install(FILES ${headers} DESTINATION include/bifrost)
install(FILES xxhash.c DESTINATION include/bifrost)
