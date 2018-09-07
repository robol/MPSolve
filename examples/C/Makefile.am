if BUILD_EXAMPLES

bin_PROGRAMS = root_of_unity

root_of_unity_SOURCES = \
	root_of_unity.c \
	dummy.cpp

root_of_unity_CFLAGS = \
	-I${top_srcdir}/include \
	-I${top_builddir}/include

root_of_unity_LDADD = \
	${top_builddir}/src/libmps/libmps.la \
	$(GMP_LIBS) 

endif
