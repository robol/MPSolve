set(TEST_CPP_FILES
  check_formal)

set(TEST_C_FILES
  check_convex check_context check_mpc
  check_matrix check_dpe 
  check_multithread check_cluster check_chebyshev check_parser check_utils 
  check_monomial_poly check_list check_secsolve check_unisolve)

if (MPSOLVE_BUILD_TESTS)
  if (CHECK_FOUND)

    add_definitions("-D_MPS_PRIVATE")
    include_directories(
      ${CMAKE_SOURCE_DIR}/include ${CMAKE_BINARY_DIR}/include
      ${GMP_INCLUDE_DIRS} ${CHECK_INCLUDE_DIRS}
      ${CMAKE_CURRENT_SOURCE_DIR})

    foreach (test_file ${TEST_CPP_FILES})
      add_executable(${test_file} "${test_file}.cpp" check_implementation.c)
      target_link_libraries(${test_file} ${GMP_LIBRARIES}
	${GMPXX_LIBRARIES} ${CHECK_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT} mps-private)
      add_test(${test_file} ${test_file})	  
    endforeach (test_file)
    
    foreach (test_file ${TEST_C_FILES})
      add_executable(${test_file} "${test_file}.c" check_implementation.c)
      target_link_libraries(${test_file} ${GMP_LIBRARIES}
	${GMPXX_LIBRARIES} ${CHECK_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT} mps-private)
      add_test(${test_file} ${test_file})
      set_tests_properties(${test_file}
	PROPERTIES ENVIRONMENT "srcdir=${CMAKE_CURRENT_SOURCE_DIR}")
    endforeach (test_file)
    

    
    
  endif (CHECK_FOUND)    
endif (MPSOLVE_BUILD_TESTS)
