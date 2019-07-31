# Same usual version stuff at the top.

cmake_minimum_required(VERSION 3.1)
project(LALIBE)

# Here I am globbing .cc stuff.

file(GLOB LIB_SRC RELATIVE ${PROJECT_SOURCE_DIR}
    "measurements/*.cc"
    "contractions/*.cc"
    "matrix_elements/*.cc"
    "io/*.cc"
    "numerics/*.cc"
    "momentum/*.cc"
)

# Build with hdf5 io stuff if the user wants.
if(BUILD_HDF5)
  add_definitions(-DBUILD_HDF5)
endif(BUILD_HDF5)

# FFTW library was found, so turn on source code that needs it.
if(FFTW_INSTALL)
  add_definitions(-DBUILD_FFTW)
endif(FFTW_INSTALL)

# Install the library: liblalibe...liblb for short.

add_library(lb ${LIB_SRC})
target_include_directories(lb INTERFACE ${CMAKE_SOURCE_DIR}/lib)
install (FILES ${CMAKE_CURRENT_BINARY_DIR}/liblb.a DESTINATION lib)

# The header files have to go somewhere. ~_~

foreach( loop_var inline hadron glue io )
  file(GLOB  H_SRC RELATIVE ${PROJECT_SOURCE_DIR} ${loop_var}/*.h )
  install (FILES ${H_SRC} DESTINATION include/${loop_var})
endforeach(loop_var)
