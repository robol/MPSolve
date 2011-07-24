if HAVE_DOXYGEN

doxyfile.stamp:
	$(DOXYGEN) Doxyfile
	echo "Timestamp" > doxyfile.stamp

CLEANFILES = doxyfile.stamp
EXTRA_DIST = Doxyfile doxygen/doxygen.css

all-local: doxyfile.stamp

clean-local:
	rm -rf $(top_builddir)/doc/html
	rm -rf $(top_builddir)/doc/latex

endif
