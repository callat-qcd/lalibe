target_include_directories( lalibe
   PRIVATE
   ${CMAKE_CURRENT_LIST_DIR})
# Sources in subdirs
add_subdirectory(measurements)
add_subdirectory(contractions)
add_subdirectory(matrix_elements)
add_subdirectory(io)
add_subdirectory(numerics)
add_subdirectory(momentum)
add_subdirectory(util)
