NULL = 

bin_PROGRAMS = secsolve

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
