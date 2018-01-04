file(GLOB sources *.cpp *.c)
file(GLOB headers *.h *.hpp *.tpp)

list(REMOVE_ITEM sources Bifrost.cpp)

add_definitions( -DMAX_KMER_SIZE=${MAX_KMER_SIZE} )

add_library(bifrost_static STATIC ${sources} ${headers})
add_library(bifrost SHARED ${sources} ${headers})

target_include_directories(bifrost_static PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_include_directories(bifrost PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

add_executable(Bifrost Bifrost.cpp)

find_package( ZLIB REQUIRED )
target_link_libraries( bifrost_static ${ZLIB_LIBRARIES})
target_link_libraries( bifrost ${ZLIB_LIBRARIES})

if ( ZLIB_FOUND )
    include_directories( ${ZLIB_INCLUDE_DIRS} )
else()
    message(FATAL_ERROR "zlib not found. Required for to output files" )
endif( ZLIB_FOUND )

find_library(JEMALLOC libjemalloc_pic.a)
if ( JEMALLOC )
	target_link_libraries( bifrost_static ${JEMALLOC} -ldl)
	target_link_libraries( bifrost ${JEMALLOC} -ldl)
else()
	target_link_libraries( bifrost_static jemalloc)
	target_link_libraries( bifrost jemalloc)
endif( JEMALLOC )

#find_library(CROARING libroaring.a)
#if ( CROARING )
#	target_link_libraries( bifrost_static ${CROARING})
#	target_link_libraries( bifrost ${CROARING})
#else()
#	target_link_libraries( bifrost_static roaring)
#	target_link_libraries( bifrost roaring)
#endif( CROARING )

find_package( Threads REQUIRED )
target_link_libraries(bifrost_static ${CMAKE_THREAD_LIBS_INIT})
target_link_libraries(bifrost ${CMAKE_THREAD_LIBS_INIT})

target_link_libraries(Bifrost bifrost)

install(TARGETS Bifrost DESTINATION bin)
install(TARGETS bifrost DESTINATION lib)
install(TARGETS bifrost_static DESTINATION lib)
install(FILES ${headers} DESTINATION include/bifrost)
