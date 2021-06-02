target_include_directories( lalibe
   PRIVATE
   ${CMAKE_CURRENT_LIST_DIR})

target_sources(lalibe
   PRIVATE
   hdf5_write_obj_funcmap.cc
   hdf5_write_obj_funcmap.h
   lalibe_qprop_io.cc
   lalibe_qprop_io.h )
