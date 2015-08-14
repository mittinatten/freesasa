if COND_GCOV
AM_CFLAGS = --coverage
AM_LDFLAGS = -lgcov
GCOV_FILES = *.gcda *.gcno *.gcov
else
GCOV_FILES =
endif

bin_PROGRAMS = freesasa
lib_LIBRARIES = libfreesasa.a
noinst_PROGRAMS = example
include_HEADERS = freesasa.h
libfreesasa_a_SOURCES = classify.c classify_user.c classify.h \
	coord.c coord.h pdb.c pdb.h \
	sasa_lr.c sasa_sr.c sasa.h \
	freesasa.c freesasa.h srp.c srp.h structure.c structure.h \
	adjacency.h adjacency.c
freesasa_SOURCES = main.c 
example_SOURCES = example.c
freesasa_LDADD = libfreesasa.a
example_LDADD = libfreesasa.a


CLEANFILES = $(GCOV_FILES)
