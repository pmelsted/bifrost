file(GLOB sources *.cpp roaring.c)
file(GLOB headers *.h *.hpp *.hh *.tcc)

list(REMOVE_ITEM sources Bifrost.cpp)

add_definitions(-DMAX_KMER_SIZE=${MAX_KMER_SIZE})

if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    # Make sure to add trailing space
    set(PROFILE_FLAGS "$<$<CONFIG:Profile>:-pg> ")
else()
    if("${CONFIG}" STREQUAL "Profile" OR "${CMAKE_BUILD_TYPE}" STREQUAL "Profile")
        message(FATAL_ERROR "Profiling is only supported with GCC")
    endif()
endif()

set(BIFROST_COMPILER_FLAGS "${PROFILE_FLAGS}-march=native")
set(BIFROST_LINKER_FLAGS "${PROFILE_FLAGS}")
set(BIFROST_COMPILER_DEFINITIONS "XXH_NAMESPACE=BIFROST_HASH_")

add_library(bifrost_static STATIC ${sources} ${headers})
add_library(bifrost_dynamic SHARED ${sources} ${headers})

set_target_properties(bifrost_static PROPERTIES OUTPUT_NAME "bifrost")
set_target_properties(bifrost_dynamic PROPERTIES OUTPUT_NAME "bifrost")

target_include_directories(bifrost_static PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_include_directories(bifrost_dynamic PUBLIC
    $<INSTALL_INTERFACE:include>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    )

target_compile_features(bifrost_static PUBLIC cxx_std_11)
target_compile_features(bifrost_static PRIVATE c_std_11)
target_compile_features(bifrost_dynamic PUBLIC cxx_std_11 c_std_11)
target_compile_features(bifrost_dynamic PRIVATE c_std_11)

target_compile_options(bifrost_static PUBLIC "${BIFROST_COMPILER_FLAGS}")
target_compile_options(bifrost_dynamic PUBLIC "${BIFROST_COMPILER_FLAGS}")

target_link_options(bifrost_static PUBLIC "${BIFROST_LINKER_FLAGS}")
target_link_options(bifrost_dynamic PUBLIC "${BIFROST_LINKER_FLAGS}")

target_compile_definitions(bifrost_static PRIVATE "${BIFROST_COMPILER_DEFINITIONS}")
target_compile_definitions(bifrost_dynamic PRIVATE "${BIFROST_COMPILER_DEFINITIONS}")

find_package(Threads REQUIRED)
target_link_libraries(bifrost_static PUBLIC ${CMAKE_THREAD_LIBS_INIT})
target_link_libraries(bifrost_dynamic PUBLIC ${CMAKE_THREAD_LIBS_INIT})

find_package(ZLIB REQUIRED)

if (ZLIB_FOUND)
    target_link_libraries(bifrost_static PUBLIC ${ZLIB_LIBRARIES})
    target_link_libraries(bifrost_dynamic PUBLIC ${ZLIB_LIBRARIES})

    target_include_directories(bifrost_static PUBLIC ${ZLIB_INCLUDE_DIRS})
    target_include_directories(bifrost_dynamic PUBLIC ${ZLIB_INCLUDE_DIRS})
else()
    message(FATAL_ERROR "zlib not found. Required for to output files")
endif(ZLIB_FOUND)

add_executable(Bifrost Bifrost.cpp)
set_property(SOURCE BlockedBloomFilter.cpp APPEND_STRING PROPERTY COMPILE_FLAGS " -funroll-loops")

target_include_directories(Bifrost PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

target_compile_features(Bifrost PUBLIC cxx_std_11 c_std_11)
target_compile_options(Bifrost PUBLIC "${BIFROST_COMPILER_FLAGS}")
target_link_options(Bifrost PUBLIC "${BIFROST_LINKER_FLAGS}")
target_compile_definitions(Bifrost PRIVATE "${BIFROST_COMPILER_DEFINITIONS}")
target_link_libraries(Bifrost PUBLIC bifrost_dynamic)

install(TARGETS Bifrost DESTINATION bin)
install(TARGETS bifrost_dynamic DESTINATION lib)
install(TARGETS bifrost_static DESTINATION lib)
install(FILES ${headers} DESTINATION include/bifrost)
install(FILES xxhash.c DESTINATION include/bifrost)
