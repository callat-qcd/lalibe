cmake_minimum_required(VERSION 3.14)
project(LALIBE
   VERSION 1.0.0.0
   LANGUAGES CXX )

# Add path to Lalibe CMake modules
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR}/configuration)


# Build Options
option(BUILD_HDF5 "build with hdf5 reading/writing" OFF)

# Find dependencies
find_package(CHROMA)

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
find_package(CHROMA)
find_package(OpenMP)
target_link_libraries(lalibe PRIVATE CHROMA)
target_link_libraries(lalibe PRIVATE OpenMP::OpenMP_CXX)

# Compile definitions
if (BUILD_HDF5)
   target_compile_definitions(lalibe PRIVATE BUILD_HDF5)
endif ()

# Install executable
install(TARGETS lalibe DESTINATION bin)
