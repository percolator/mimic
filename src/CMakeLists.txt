
#cmake -DCMAKE_INSTALL_PREFIX=/tmp/install /home/lukask/cvs/percolator/mimic/src && make && make install
cmake_minimum_required(VERSION 2.6)



# configure_file(${CMAKE_CURRENT_SOURCE_DIR}/config.h.cmake ${CMAKE_CURRENT_BINARY_DIR}/config.h)


include_directories(${CMAKE_CURRENT_SOURCE_DIR})


# set(R_SRCS stringa.F stringch.F numrec.F anal_main.F probfcn.F vector.F vdw.F alignanal.F io.F profiles.F )



# add_library(trace STATIC ${TRACE_SRCS} )

add_executable(mimic AminoAcidDist.cpp  Option.cpp Peptides.cpp)
target_link_libraries(mimic m)
install(TARGETS mimic DESTINATION bin)


#set_target_properties(mimic PROPERTIES COMPILE_FLAGS "${extraflag} ")
