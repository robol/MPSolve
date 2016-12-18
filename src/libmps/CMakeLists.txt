#
# Testing the basics
#
include(CheckFunctionExists)
include(CheckCSourceCompiles)

check_function_exists(getline HAVE_GETLINE)

# Check for the visibility attribute
check_c_source_compiles(
  "int __attribute__ ((visibility (\"hidden\"))) main (void) { return 1; }"
  HAVE_HIDDEN_VISIBILITY_ATTRIBUTE)

# Check for support of complex numbers
check_c_source_compiles(
  "#include <complex.h> 
  int main() { 
    double complex a, b;
  }"
  HAVE_COMPLEX_H)

if (NOT HAVE_COMPLEX_H)
  set (MPS_USE_BUILTIN_COMPLEX 0 CACHE INTERNAL MPS_USE_BUILTIN_COMPLEX)
else(NOT HAVE_COMPLEX_H)
  set (MPS_USE_BUILTIN_COMPLEX 1 CACHE INTERNAL MPS_USE_BUILTIN_COMPLEX)
endif (NOT HAVE_COMPLEX_H)

#
# Source files for libmps
#
set(LIBMPS_SOURCES
	chebyshev/chebyshev-evaluation.c 
	chebyshev/chebyshev-parser.c 
	chebyshev/chebyshev.c 
	common/aberth.c 
	common/approximation.c 
	common/cluster-analysis.c 
	common/cluster.c 
	common/context.c 
	common/convex.c 
	common/defaults.c 
	common/file-starting.c 
	common/improve.c 
	common/inclusion.c 
	common/inline-poly-parser.c 
	common/interface.c 
	common/jacobi-aberth.c 
	common/list-element.c 
	common/list.c 
	common/modify.c 
	common/nroots-polynomial.cpp 
	common/parser.c 
	common/polynomial.c 
	common/polynomialxx.cpp 
	common/recursive-starting.c 
	common/sort.c 
	common/starting-configuration.c 
	common/starting.c 
        common/strndup.c 
	common/test.c 
	common/tools.c 
	common/touch.c 
	common/user.c 
	common/utils.c 
	common/validation.c 
	formal/formal-monomial.cpp 
	formal/formal-polynomial.cpp 
	floating-point/gmptools.c 
	floating-point/link.c 
	floating-point/mpc.c 
	floating-point/mt.c 
	general/general-radius.c 
	general/general-starting.c 
	matrix/hessenberg-determinant.c 
	monomial/horner.c 
	monomial/monomial-matrix-poly.c 
	monomial/monomial-parser.c 
	monomial/monomial-poly.c 
	monomial/monomial-threading.c 
	monomial/newton.c 
	${CMAKE_BINARY_DIR}/src/libmps/yacc-parser.c
	${CMAKE_BINARY_DIR}/src/libmps/tokenizer.c
	monomial/shift.c 
	secsolve/secular-ga.c 
	secsolve/secular-iteration.c 
	secsolve/secular-regeneration.c 
	secsolve/standard-regeneration-driver.c 
	secular/secular-equation.c 
	secular/secular-evaluation.c 
	secular/secular-newton.c 
	secular/secular-parser.c 
	secular/secular-starting.c 
	system/abstract-input-stream.cpp 
	system/file-input-stream.cpp 
	system/memory-file-stream.cpp 
	system/data.c 
	system/debug.c 
	system/getline.c 
	system/getopts.c 
	system/input-buffer.c 
	system/input-output.c 
	system/threading.c 
	unisolve/main.c
	unisolve/solve.c)

find_program(BISON NAMES bison)
find_program(FLEX NAMES flex)

add_custom_command(
  COMMAND ${FLEX}
  ARGS -o ${CMAKE_BINARY_DIR}/src/libmps/tokenizer.c
       ${CMAKE_SOURCE_DIR}/src/libmps/monomial/tokenizer.l
  OUTPUT ${CMAKE_BINARY_DIR}/src/libmps/tokenizer.c)

add_custom_command(
  COMMAND ${BISON}
  ARGS -y ${CMAKE_SOURCE_DIR}/src/libmps/monomial/yacc-parser.y
       -o ${CMAKE_BINARY_DIR}/src/libmps/yacc-parser.c
       -d 
  DEPENDS ${CMAKE_BINARY_DIR}/src/libmps/tokenizer.c
  OUTPUT ${CMAKE_BINARY_DIR}/src/libmps/yacc-parser.c)

set_source_files_properties(
  ${CMAKE_BINARY_DIR}/src/libmps/monomial/yacc-parser.c GENERATED)
set_source_files_properties(
  ${CMAKE_BINARY_DIR}/src/libmps/monomial/tokenizer.c GENERATED)

#
# Checking for mandatory libraries required to build libmps
#
find_package(GMP REQUIRED)

#
# Checking for optional libraries required to build libmps
#
find_package(Tcmalloc)
find_package(PkgConfig)
if (PkgConfig_FOUND)
  
  pkg_check_modules(GTK3 gtk+-3.0)
  if (GTK3_FOUND)
    # This is needed for the config.h file
    set(HAVE_GTK 1 CACHE INTERNAL HAVE_GTK)
  endif (GTK3_FOUND)
  
endif (PkgConfig_FOUND)

add_library(mps SHARED
	${LIBMPS_SOURCES})

target_link_libraries(mps
  ${GMP_LIBRARIES} ${GMPXX_LIBRARIES}
  ${TCMALLOC_LIBRARY} ${GTK3_LIBRARIES})

link_directories(${GTK3_LIBRARY_DIRS})

add_definitions("-D_MPS_PRIVATE -include ${CMAKE_BINARY_DIR}/config.h")

configure_file(${CMAKE_SOURCE_DIR}/include/mps/mt.h.in
  ${CMAKE_BINARY_DIR}/include/mps/mt.h)

target_include_directories(mps
  PUBLIC ${CMAKE_SOURCE_DIR}/include ${CMAKE_BINARY_DIR}/include
  PRIVATE ${GMP_INCLUDE_DIR} ${TCMALLOC_INCLUDE_DIR} ${GTK3_INCLUDE_DIRS})

