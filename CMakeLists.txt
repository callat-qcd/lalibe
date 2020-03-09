# Versions of cmake and code set here.

cmake_minimum_required(VERSION 3.1)
project(LALIBE)
set (LALIBE_VERSION_MAJOR 1)
set (LALIBE_VERSION_MINOR 0)
set (LALIBE_VERSION_PATCH 0)
set (LALIBE_VERSION_TWEAK 0)

# Usual cmake header configure stuff.

configure_file (
  "${PROJECT_SOURCE_DIR}/configuration/lalibe_configure.h.in"
  "${PROJECT_BINARY_DIR}/lalibe_configure.h"
  )

# Have an on/off switch for hdf5.
set(BUILD_HDF5 OFF CACHE BOOL "build with hdf5 reading/writing")

# Append cmake stuff with chroma installation path.

list(APPEND CMAKE_PREFIX_PATH ${CHROMA_INSTALL})
message("I hope for your sake that chroma is installed here: " ${CHROMA_INSTALL})

# Find chroma configuration info so we can borrow it for our own install.

find_package(CHROMA PATHS configuration)
set(CMAKE_CXX_FLAGS "${CHROMA_CXX_FLAGS} -I${CMAKE_SOURCE_DIR}/lib")
add_subdirectory(lib)
set(ORDERED_LIBS "-L${CMAKE_CURRENT_BINARY_DIR}/lib -llb ${CHROMA_LD_FLAGS} ${CHROMA_LIBS}")
message("Here are the chroma cxx flags I found: " ${CHROMA_CXX_FLAGS})
message("Here are the chroma ld flags I found: " ${CHROMA_LD_FLAGS})
message("Here are the chroma libs I found: " ${CHROMA_LIBS})

#These lines are usually not needed, but may be useful for debugging purposes.
find_package(MPI REQUIRED)
include_directories(${MPI_INCLUDE_PATH})
message("Here is where I think mpi is:  " ${MPI_INCLUDE_PATH})
#set(CMAKE_CXX_COMPILE_FLAGS ${CMAKE_CXX_COMPILE_FLAGS} ${MPI_COMPILE_FLAGS})
#set(CMAKE_CXX_LINK_FLAGS ${CMAKE_CXX_LINK_FLAGS} ${MPI_LINK_FLAGS})
#message("Here is what I think the mpi compile flags should be: " ${MPI_COMPILE_FLAGS})
#message("Here is what I think the mpi is link flags should be: " ${MPI_LINK_FLAGS})

# Time for the magic to happen.

add_executable(lalibe main/lalibe.cc)
add_dependencies(lalibe lb)
target_include_directories(lalibe INTERFACE ${CMAKE_SOURCE_DIR}/lib)

# Zelda can only be jealous of how our files Link.

target_link_libraries(lalibe ${ORDERED_LIBS}) 
install (TARGETS lalibe DESTINATION bin)



