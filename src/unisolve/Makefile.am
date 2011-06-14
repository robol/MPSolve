NULL = 

bin_PROGRAMS = unisolve

unisolve_CFLAGS = \
	$(UNISOLVE_CFLAGS) \
	$(NULL)

unisolve_SOURCES = \
	unisolve.c \
	$(NULL)

unisolve_LDADD = \
	-L../libmps \
	-lmps \
	-lgmp \
	-lm \
	$(NULL)
