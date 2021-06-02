target_include_directories( lalibe
   PRIVATE
   ${CMAKE_CURRENT_LIST_DIR})

target_sources(lalibe
   PRIVATE
    spin_basis.h
    spin_basis.cc)
