target_include_directories( lalibe
   PRIVATE
   ${CMAKE_CURRENT_LIST_DIR})

target_sources(lalibe
   PRIVATE
   baryon_contractions_func_w.cc
   baryon_seqsource_w.h
   meson_contractions_func_w.h
   seqsource_contractions_func_w.cc
   baryon_contractions_func_w.h
   proton_contractions_func_w.cc
   seqsource_contractions_func_w.h
   baryon_seqsource_w.cc
   meson_contractions_func_w.cc
   proton_contractions_func_w.h 
    pipi_scattering_func_w.cc
    pipi_scattering_func_w.h
)

