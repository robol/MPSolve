TESTS = simple_polynomial.py algorithms.py

# Make sure that we can load the python module from the test environment and that
# the locally built module is always preferred to the system one. 
AM_TESTS_ENVIRONMENT = \
	export LD_LIBRARY_PATH=${top_builddir}/src/libmps/.libs; \
	export PYTHONPATH=${top_srcdir}/examples/python/;
	
EXTRA_DIST = simple_polynomial.py algorithms.py
