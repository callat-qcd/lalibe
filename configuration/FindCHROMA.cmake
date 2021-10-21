#[=======================================================================[:
FindCHROMA
-------

Finds the CHROMA library.

Imported Targets
^^^^^^^^^^^^^^^^

This module provides the following imported target, if found:

``CHROMA``
  The CHROMA library

Result Variables
^^^^^^^^^^^^^^^^

This will define the following variables:

``CHROMA_FOUND``
  True if the CHROMA library has been found.
``CHROMA_CXX_FLAGS``
  The CXX flags of the CHROMA library which was found.
``CHROMA_INCLUDE_DIRS``
  Include directories needed to use CHROMA.
``CHROMA_LD_FLAGS``
  Linker flags needed to link to CHROMA.
``CHROMA_LIBRARIES``
  Libraries needed to link to CHROMA.
#]=======================================================================]

find_program( _chroma_config_exe
   NAMES chroma-config
   PATHS ${CHROMA_ROOT}/bin ENV CHROMA_ROOT
   )

if ( _chroma_config_exe )
   set(CHROMA_INCLUDE_DIRS ${CHROMA_ROOT}/include)
   execute_process(COMMAND ${_chroma_config_exe} --version
      OUTPUT_VARIABLE CHROMA_VERSION
      OUTPUT_STRIP_TRAILING_WHITESPACE)
   execute_process(COMMAND ${_chroma_config_exe} --cxxflags
      OUTPUT_VARIABLE CHROMA_CXX_FLAGS
      OUTPUT_STRIP_TRAILING_WHITESPACE)
   execute_process(COMMAND ${_chroma_config_exe} --ldflags
      OUTPUT_VARIABLE CHROMA_LD_FLAGS
      OUTPUT_STRIP_TRAILING_WHITESPACE)
   execute_process(COMMAND ${_chroma_config_exe} --libs
      OUTPUT_VARIABLE CHROMA_LIBRARIES
      OUTPUT_STRIP_TRAILING_WHITESPACE)
endif ()

include(FindPackageHandleStandardArgs)

find_package_handle_standard_args(CHROMA
   REQUIRED_VARS
   CHROMA_LIBRARIES
   CHROMA_INCLUDE_DIRS
   CHROMA_LD_FLAGS
   CHROMA_CXX_FLAGS
   VERSION_VAR
   CHROMA_VERSION
   )

mark_as_advanced(CHROMA_VERSION CHROMA_LIBRARIES CHROMA_CXX_FLAGS
   CHROMA_INCLUDE_DIRS CHROMA_LD_FLAGS)

# Create imported target
if (CHROMA_FOUND AND NOT TARGET CHROMA)
    if(CHROMA_CXX_FLAGS)
        string(REPLACE " " ";" CHROMA_CXX_FLAGS "${CHROMA_CXX_FLAGS}")
       endif()
     if(CHROMA_LD_FLAGS)
        string(REPLACE " " ";" CHROMA_LD_FLAGS "${CHROMA_LD_FLAGS}")
    endif()
    # add_library(CHROMA UNKNOWN IMPORTED GLOBAL)
    add_library(CHROMA INTERFACE)

    # If we build QUDA and CHROMA with QIO, then Lalibe has difficulty 
    # linking with QIO/LIME unless we re-specify the linking
    if(QIO_LINK)
        set_target_properties(CHROMA
          PROPERTIES
          INTERFACE_COMPILE_OPTIONS "${CHROMA_CXX_FLAGS}"
          INTERFACE_INCLUDE_DIRECTORIES "${CHROMA_INCLUDE_DIRS}"
          INTERFACE_LINK_OPTIONS "${CHROMA_LD_FLAGS}"
          INTERFACE_LINK_LIBRARIES "${CHROMA_LIBRARIES} -lqio -llime"
          )
    else()
        set_target_properties(CHROMA
          PROPERTIES
          INTERFACE_COMPILE_OPTIONS "${CHROMA_CXX_FLAGS}"
          INTERFACE_INCLUDE_DIRECTORIES "${CHROMA_INCLUDE_DIRS}"
          INTERFACE_LINK_OPTIONS "${CHROMA_LD_FLAGS}"
          INTERFACE_LINK_LIBRARIES "${CHROMA_LIBRARIES}"
          )
    endif()
endif ()
