cmake_minimum_required(VERSION 3.0)

# This is the main configuration file for the MPSolve project.
project(MPSolve C CXX)

set(MPS_VERSION_MAJOR 3)
set(MPS_VERSION_MINOR 1)
set(MPS_VERSION_PATCH 5)

set(MPS_VERSION_STRING
  ${MPS_VERSION_MAJOR}.${MPS_VERSION_MINOR}.${MPS_VERSION_PATCH})
set(PACKAGE_STRING
  "MPSolve ${MPS_VERSION_STRING}" CACHE INTERNAL PACKAGE_STRING)

enable_testing()

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH}
  "${CMAKE_SOURCE_DIR}/cmake/")

# Project-wide includes 
include_directories(
  ${CMAKE_SOURCE_DIR}/include
  ${CMAKE_BINARY_DIR}/include)

add_subdirectory(src)
add_subdirectory(examples)

# Create the config.h file that was needed by the old automake setup.
# As of now, we still use the file to track many variables. 
configure_file(${CMAKE_SOURCE_DIR}/cmake/config.h.in
  ${CMAKE_BINARY_DIR}/include/mps/config.h)
