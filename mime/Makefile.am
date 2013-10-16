# Inspired by Homebank MIME installation Makefile.am

mimedir = $(datadir)/mime-info
mime_DATA = mpsolve.mime

mimexmldir = $(datadir)/mime/packages
mimexml_DATA = mpsolve.xml

EXTRA_DIST = \
	$(mime_DATA) \
	$(mimexml_DATA) 

## Update mime only when system wide installation takes place

install-data-local: install-mimexmlDATA
	if [ -f $(DESTDIR)$(datadir)/mime/packages/freedesktop.org.xml ] ; then \
	  if which update-mime-database>/dev/null 2>&1; then \
	    update-mime-database $(DESTDIR)$(datadir)/mime; \
	  fi \
	fi

uninstall-local:
	if [ -f $(DESTDIR)$(datadir)/mime/packages/freedesktop.org.xml ] ; then \
	  if which update-mime-database>/dev/null 2>&1; then \
	    update-mime-database $(DESTDIR)$(datadir)/mime; \
	  fi \
	fi
