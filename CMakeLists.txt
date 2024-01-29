cmake_minimum_required(VERSION 3.14)
project(LALIBE
   VERSION 1.0.0.0
   LANGUAGES CXX C)

# Add path to Lalibe CMake modules
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR}/configuration)


# Build Options
option(BUILD_HDF5 "build with hdf5 reading/writing" OFF)

set(Chroma_DIR ${CHROMA_ROOT})
string(FIND ${CHROMA_ROOT} "cmake" CMAKE_BUILD)

# Find dependencies
if( ${CMAKE_BUILD} EQUAL "-1" )
  message("Building with autotools\n")
  find_package(CHROMA)
else()
  find_package(Chroma REQUIRED)
endif()

# Create lalibe target and add to it main source and
# config header
add_executable(lalibe)

configure_file (
  "${PROJECT_SOURCE_DIR}/configuration/lalibe_configure.h.in"
  "${PROJECT_BINARY_DIR}/lalibe_configure.h"
  )

target_sources(lalibe
   PRIVATE
   main/lalibe.cc
   ${PROJECT_BINARY_DIR}/lalibe_configure.h
   )

target_include_directories(lalibe PRIVATE ${PROJECT_BINARY_DIR})

# Add all the remaining sources
add_subdirectory(lib)

# Find dependencies and link them to lalibe
if( ${CMAKE_BUILD} EQUAL "-1" )
  find_package(CHROMA)
  if ( CHROMA_FOUND )
      message("\nCHROMA_ROOT = ${CHROMA_ROOT}\n")
  else()
      message(FATAL_ERROR "we could not find CHROMA_ROOT: ${CHROMA_ROOT}")
  endif()
  target_link_libraries(lalibe PRIVATE CHROMA)
  target_link_libraries(lalibe PRIVATE qio lime)
else()
  target_link_libraries(lalibe PRIVATE Chroma::chromalib)
endif()

find_package(OpenMP)
target_link_libraries(lalibe PRIVATE OpenMP::OpenMP_CXX)

# Compile definitions
if (BUILD_HDF5)
   target_compile_definitions(lalibe PRIVATE BUILD_HDF5)
endif ()

# Install executable
install(TARGETS lalibe DESTINATION bin)
