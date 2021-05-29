target_include_directories( lalibe
   PRIVATE
   ${CMAKE_CURRENT_LIST_DIR})

target_sources(lalibe
   PRIVATE
   binaryRecursiveColoring.cc
   binaryRecursiveColoring.h
   binaryRecursiveColoring_v2.h)
