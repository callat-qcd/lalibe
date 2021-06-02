target_include_directories( lalibe
   PRIVATE
   ${CMAKE_CURRENT_LIST_DIR})

target_sources(lalibe
   PRIVATE
   bilinear_gamma.cc
   chromomag_seqsource_w.cc
   lalibe_formfac_w.h
   tr_less_fields.cc
   bilinear_gamma.h
   chromomag_seqsource_w.h
   lalibe_formfac_w.cc
   tr_less_fields.h )
