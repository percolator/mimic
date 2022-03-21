project(mimic)
#cmake -DCMAKE_INSTALL_PREFIX=/tmp/install /home/lukask/cvs/percolator/mimic/src && make && make install
cmake_minimum_required(VERSION 3.18)

# configure_file(${CMAKE_CURRENT_SOURCE_DIR}/config.h.cmake ${CMAKE_CURRENT_BINARY_DIR}/config.h)


include_directories(${CMAKE_CURRENT_SOURCE_DIR})


# set(R_SRCS stringa.F stringch.F numrec.F anal_main.F probfcn.F vector.F vdw.F alignanal.F io.F profiles.F )



# add_library(trace STATIC ${TRACE_SRCS} )

add_executable(mimic mimic.cpp AminoAcidDist.cpp  Option.cpp Peptides.cpp)
target_link_libraries(mimic m)
target_compile_features(mimic PRIVATE cxx_std_17)
install(TARGETS mimic DESTINATION bin)

#set_target_properties(mimic PROPERTIES COMPILE_FLAGS "${extraflag} ")

# Download and run conan.cmake
if (NOT EXISTS "${CMAKE_BINARY_DIR}/conan.cmake")
    message(STATUS "Downloading conan.cmake from https://github.com/conan-io/cmake-conan")
    file(DOWNLOAD "https://raw.githubusercontent.com/conan-io/cmake-conan/v0.16.1/conan.cmake"
            "${CMAKE_BINARY_DIR}/conan.cmake"
            EXPECTED_HASH SHA256=396e16d0f5eabdc6a14afddbcfff62a54a7ee75c6da23f32f7a31bc85db23484
            TLS_VERIFY ON)
endif ()
include(${CMAKE_BINARY_DIR}/conan.cmake)
conan_cmake_run(CONANFILE conanfile_gcc.txt
        PROFILE gcc
        BASIC_SETUP CMAKE_TARGETS
        BUILD missing)

include(${CMAKE_BINARY_DIR}/conanbuildinfo.cmake)
# These tests can use the Catch2-provided main
add_executable(tests test.cpp AminoAcidDist.cpp  Option.cpp Peptides.cpp)
target_compile_features(tests PRIVATE cxx_std_17)
target_link_libraries(tests PRIVATE CONAN_PKG::catch2)
configure_file(${PROJECT_SOURCE_DIR}/../min.fasta ${CMAKE_BINARY_DIR}/bin/ COPYONLY)
