NULL = 

bin_PROGRAMS = secsolve
check_PROGRAMS = secsolve-test

secsolve_CFLAGS = \
	$(SECSOLVE_CFLAGS) \
	$(NULL)

secsolve_SOURCES = \
	secsolve.c \
	$(NULL)

secsolve_LDADD = \
	-L../libmps \
	-lmps \
	-lgmp \
	-lm \
	$(NULL)

secsolve_test_CFLAGS = \
	$(SECSOLVE_CFLAGS) \
	$(NULL)

secsolve_test_SOURCES = \
	secsolve-test.c \
	$(NULL)

secsolve_test_LDADD = \
	-L../libmps \
	-lmps \
	-lgmp \
	-lm \
	$(NULL)
