if BUILD_EXAMPLES

bin_PROGRAMS = root_of_unity secular

root_of_unity_SOURCES = \
	root_of_unity.c \
	dummy.cpp

root_of_unity_CFLAGS = \
	-I${top_srcdir}/include \
	-I${top_builddir}/include

root_of_unity_LDADD = \
	${top_builddir}/src/libmps/libmps.la \
	$(GMP_LIBS) 


secular_SOURCES = \
	secular.c \
	dummy.cpp

secular_CFLAGS = \
	-I${top_srcdir}/include \
	-I${top_builddir}/include

secular_LDADD = \
	${top_builddir}/src/libmps/libmps.la \
	$(GMP_LIBS) 

endif
