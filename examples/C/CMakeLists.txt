include_directories(
  ${CMAKE_SOURCE_DIR}/include
  ${CMAKE_BINARY_DIR}/include)

link_libraries(mps gmp m)

add_executable(root_of_unity
  root_of_unity.c)
