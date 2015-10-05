if COND_GCOV
AM_CFLAGS = --coverage
AM_LDFLAGS = -lgcov
endif
GCOV_FILES = *.gcda *.gcno *.gcov

bin_PROGRAMS = freesasa
lib_LIBRARIES = libfreesasa.a
noinst_PROGRAMS = example
include_HEADERS = freesasa.h
libfreesasa_a_SOURCES = classify.c user_config.c classify.h \
	coord.c coord.h pdb.c pdb.h \
	sasa_lr.c sasa_sr.c sasa.h \
	freesasa.c freesasa.h structure.c structure.h \
	nb.h nb.c util.h util.c \
	selector.h selector.c $(lp_output)
freesasa_SOURCES = main.c 
example_SOURCES = example.c
freesasa_LDADD = libfreesasa.a
example_LDADD = libfreesasa.a

lp_output = lexer.c lexer.h parser.c parser.h

$(lp_output): lexer.l parser.y
	flex lexer.l
	bison parser.y

CLEANFILES = $(GCOV_FILES) *~
