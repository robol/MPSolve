# Option parsing
option(MPSOLVE_DEBUG "Debug mode" ON)
option(MPSOLVE_DEBUG_BUILD "Build MPSolve without optimizations" OFF)
option(MPSOLVE_EXAMPLES "Build examples and bindings" ON)
option(MPSOLVE_GUI "Build the user interface of MPSolve" ON)
option(MPSOLVE_GRAPHICAL_DEBUGGER "Build the graphical debugger for MPSolve" ON)
option(MPSOLVE_BUILD_TESTS "Build the unit tests for MPSolve" ON)

if (MPSOLVE_BUILD_TESTS)  
  find_package(Check)
endif (MPSOLVE_BUILD_TESTS)

set(CMAKE_THREAD_PREFER_PTHREAD ON)
find_package(Threads REQUIRED)
if (NOT CMAKE_USE_PTHREADS_INIT)
  # Trigger an error
endif (NOT CMAKE_USE_PTHREADS_INIT)


add_subdirectory(libmps)
add_subdirectory(mpsolve)
add_subdirectory(xmpsolve)
add_subdirectory(tests)
