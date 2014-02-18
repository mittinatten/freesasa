if COND_GCOV
AM_CFLAGS = --coverage
AM_LDFLAGS = -lgcov
GCOV_FILES = *.gcda *.gcno *.gcov
else
GCOV_FILES =
endif

bin_PROGRAMS = calc_sasa
lib_LIBRARIES = libsasa.a
noinst_PROGRAMS = example
include_HEADERS = sasalib.h
libsasa_a_SOURCES = classify.c classify.h coord.c coord.h pdb.c pdb.h \
	sasa_lr.c sasa_sr.c sasa.h \
	sasalib.c sasalib.h srp.c srp.h structure.c structure.h
calc_sasa_SOURCES = main.c 
calc_sasa_CFLAGS = -DPROGRAM_NAME=\"calc_sasa\" $(AM_CFLAGS)
example_SOURCES = example.c 
calc_sasa_LDADD = libsasa.a
example_LDADD = libsasa.a

CLEANFILES = $(GCOV_FILES)