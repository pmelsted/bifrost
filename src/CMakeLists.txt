file(GLOB sources *.cpp)
file(GLOB headers *.h *.hpp)

list(REMOVE_ITEM sources BFGraph.cpp)

add_definitions( -DMAX_KMER_SIZE=${MAX_KMER_SIZE} )

add_library(BFGraph_core ${sources} ${headers})
target_include_directories(BFGraph_core PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

add_executable(BFGraph BFGraph.cpp)

find_package( Threads REQUIRED )
target_link_libraries(BFGraph BFGraph_core pthread)

if(LINK MATCHES static)
    set(BUILD_SHARED_LIBS OFF)

    if (UNIX AND NOT APPLE)
        #set(CMAKE_EXE_LINKER_FLAGS  "-static -static-libgcc -static-libstdc++")
        set(CMAKE_EXE_LINKER_FLAGS  "-static -static-libstdc++")
        SET(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
        set(CMAKE_EXE_LINKER_FLAGS  "-static -static-libgcc -static-libstdc++")
    endif(UNIX AND NOT APPLE)

    SET_TARGET_PROPERTIES(BFGraph PROPERTIES LINK_SEARCH_END_STATIC 1)
endif(LINK MATCHES static)

find_package( ZLIB REQUIRED )
target_link_libraries( BFGraph ${ZLIB_LIBRARIES})

if ( ZLIB_FOUND )
    include_directories( ${ZLIB_INCLUDE_DIRS} )
else()
    message(FATAL_ERROR "zlib not found. Required for to output files" )
endif( ZLIB_FOUND )



if(LINK MATCHES static)
    if (UNIX AND NOT APPLE)
        target_link_libraries(BFGraph librt.a)
    endif()
else()
    if (UNIX AND NOT APPLE)
        target_link_libraries(BFGraph rt)
    endif()
endif(LINK MATCHES static)



install(TARGETS BFGraph DESTINATION bin)