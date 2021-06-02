target_include_directories( lalibe
   PRIVATE
   ${CMAKE_CURRENT_LIST_DIR})

target_sources(lalibe
   PRIVATE
   lalibe_sftmom.cc
   lalibe_sftmom.h
   )
