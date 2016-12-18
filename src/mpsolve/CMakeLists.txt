add_executable(mpsolve
  mpsolve.c
  iteration-logger.c
  iteration-logger.h)

add_definitions("-pthread")

target_link_libraries(mpsolve
  mps pthread)
