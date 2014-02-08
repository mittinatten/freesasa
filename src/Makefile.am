AM_CFLAGS = -O2 -std=c99 -Wall -pedantic

bin_PROGRAMS = calc_sasa
lib_LIBRARIES = libsasa.a
noinst_PROGRAMS = example
include_HEADERS = sasalib.h
libsasa_a_SOURCES = classify.c classify.h coord.c coord.h pdb.c pdb.h sasa.c sasa.h \
	sasalib.c sasalib.h srp.c srp.h structure.c structure.h
calc_sasa_SOURCES = main.c 
calc_sasa_CFLAGS = -DPROGRAM_NAME=\"calc_sasa\" $(AM_CFLAGS)
example_SOURCES = example.c 
calc_sasa_LDADD = libsasa.a
example_LDADD = libsasa.a