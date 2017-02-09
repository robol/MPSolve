# Testing the basics
include(CheckFunctionExists)
include(CheckCSourceCompiles)
include(CheckCCompilerFlag)

include_directories(
  PUBLIC ${CMAKE_SOURCE_DIR}/include ${CMAKE_BINARY_DIR}/include
  PRIVATE ${GMP_INCLUDE_DIR} ${TCMALLOC_INCLUDE_DIR} ${GTK3_INCLUDE_DIRS})

check_function_exists(getline HAVE_GETLINE)
check_function_exists(strndup HAVE_STRNDUP)

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

check_c_compiler_flag("-fcx-limited-range" FLAG_CX_LIMITED_RANGE)
check_c_compiler_flag("-fno-math-errno" FLAG_NO_MATH_ERRNO)
check_c_compiler_flag("-fomit-frame-pointer" FLAG_OMIT_FRAME_POINTER)

if (FLAG_CX_LIMITED_RANGE)
  set(MPS_COMMON_DEFINITIONS "${MPS_COMMON_DEFINITIONS} -fcx-limited-range")
endif (FLAG_CX_LIMITED_RANGE)

if (FLAG_MATH_ERRNO)
  set(MPS_COMMON_DEFINITIONS "${MPS_COMMON_DEFINITIONS} -fmath-errno")
endif (FLAG_MATH_ERRNO)

if (FLAG_OMIT_FRAME_POINTER)
  set(MPS_COMMON_DEFINITIONS "${MPS_COMMON_DEFINITIONS} -fomit-frame-pointer")
endif (FLAG_OMIT_FRAME_POINTER)

# Select the right optimization flags
if (NOT MPSOLVE_DEBUG_BUILD)
  set(MPS_COMMON_DEFINITIONS "${MPS_COMMON_DEFINITIONS} -O2")
else (NOT MPSOLVE_DEBUG_BUILD)
  set(MPS_COMMON_DEFINITIONS "${MPS_COMMON_DEFINITIONS} -O0 -g")
endif (NOT MPSOLVE_DEBUG_BUILD)

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

find_package(FLEX REQUIRED)
find_package(BISON REQUIRED)

bison_target(libmps_parser monomial/yacc-parser.y
    ${CMAKE_CURRENT_BINARY_DIR}/yacc-parser.c)
flex_target(libmps_tokenizer monomial/tokenizer.l
    ${CMAKE_BINARY_DIR}/src/libmps/tokenizer.c 
    COMPILE_FLAGS "-d")
add_flex_bison_dependency(libmps_tokenizer libmps_parser)

# Checking for mandatory libraries required to build libmps
find_package(GMP REQUIRED)

# Checking for optional libraries required to build libmps
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
	${LIBMPS_SOURCES}
        ${BISON_libmps_parser_OUTPUTS}
        ${FLEX_libmps_tokenizer_OUTPUTS})

set_target_properties(mps PROPERTIES
        VERSION ${MPS_VERSION_STRING}
        SOVERSION ${MPS_VERSION_MAJOR})

target_link_libraries(mps
  ${GMP_LIBRARIES} ${GMPXX_LIBRARIES}
  ${TCMALLOC_LIBRARY} ${GTK3_LIBRARIES}
  ${CMAKE_THREAD_LIBS_INIT})

link_directories(${GTK3_LIBRARY_DIRS})

if (NOT MPSOLVE_DEBUG)
  set(MPS_COMMON_DEFINITIONS "${MPS_COMMON_DEFINITIONS} -DDISABLE_DEBUG")
else (NOT MPSOLVE_DEBUG)
  set(MPS_COMMON_DEFINITIONS "${MPS_COMMON_DEFINITIONS} -DNICE_DEBUG")
endif (NOT MPSOLVE_DEBUG)

add_definitions("${MPS_COMMON_DEFINITIONS}
   -D_MPS_PRIVATE 
   -include ${CMAKE_BINARY_DIR}/config.h")

configure_file(${CMAKE_SOURCE_DIR}/include/mps/mt.h.in
  ${CMAKE_BINARY_DIR}/include/mps/mt.h)

# Handling tests
if (CHECK_FOUND)
  # If the check library is found we build a custom version of libmps
  # that exposes all the private methods, and allow for fine-grained
  # testing.   
  add_library(mps-private SHARED
    ${LIBMPS_SOURCES}
    ${BISON_libmps_parser_OUTPUTS}
    ${FLEX_libmps_tokenizer_OUTPUTS})

  target_link_libraries(mps-private
    ${GMP_LIBRARIES} ${GMPXX_LIBRARIES}
    ${TCMALLOC_LIBRARY} ${GTK3_LIBRARIES}
    ${CMAKE_THREAD_LIBS_INIT})

  # Make sure that all methods are exposed in mps-private, since this is used
  # for testing and needs to access all the code base. 
  target_compile_definitions(mps-private
    PRIVATE "-DMPS_PUBLISH_PRIVATE_METHODS")
endif (CHECK_FOUND)

# Handling installation
install(TARGETS mps
  LIBRARY DESTINATION lib)

install(DIRECTORY ${CMAKE_SOURCE_DIR}/include/mps
  DESTINATION include
  FILES_MATCHING PATTERN "*.h")
