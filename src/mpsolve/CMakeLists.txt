include_directories(
  ${CMAKE_SOURCE_DIR}/include
  ${CMAKE_BINARY_DIR}/include)

add_executable(mpsolve
  mpsolve.c
  iteration-logger.c
  iteration-logger.h)

target_link_libraries(mpsolve
  mps ${CMAKE_THREAD_LIBS_INIT})

install(TARGETS mpsolve
  RUNTIME DESTINATION bin)
