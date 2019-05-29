# Take relevant flags from chroma configuration (call chroma-config). If we can't find chroma, we complain.

find_program(CHROMA_CONFIG NAMES chroma-config)

if( NOT CHROMA_CONFIG ) 
  message(FATAL_ERROR "Dude, I can't find where your chroma is. I am sorry, I have to exit now. :(")  
endif( NOT CHROMA_CONFIG)

execute_process(COMMAND ${CHROMA_CONFIG} --cxxflags OUTPUT_VARIABLE CHROMA_CXX_FLAGS  OUTPUT_STRIP_TRAILING_WHITESPACE)
execute_process(COMMAND ${CHROMA_CONFIG} --ldflags OUTPUT_VARIABLE CHROMA_LD_FLAGS OUTPUT_STRIP_TRAILING_WHITESPACE)
execute_process(COMMAND ${CHROMA_CONFIG} --libs  OUTPUT_VARIABLE CHROMA_LIBS OUTPUT_STRIP_TRAILING_WHITESPACE)

