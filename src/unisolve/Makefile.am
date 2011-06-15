NULL = 

bin_PROGRAMS = unisolve
check_PROGRAMS = unisolve-test

unisolve_test_SOURCES = \
	unisolve-test.c \
	$(NULL)

unisolve_test_CFLAGS = \
	$(UNISOLVE_CFLAGS) \
	$(NULL)

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

unisolve_test_LDADD = \
	-L../libmps \
	-lmps \
	-lgmp \
	-lm \
	$(NULL)
