cmake_minimum_required(VERSION 2.8)

#
# This is the main configuration file for the MPSolve project.
#
project(MPSolve)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH}
  "${CMAKE_SOURCE_DIR}/cmake/")

add_subdirectory(src)

# Create the config.h file that was needed by the old automake setup.
configure_file(${CMAKE_SOURCE_DIR}/config.h.in ${CMAKE_BINARY_DIR}/config.h)

