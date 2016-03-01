AM_LDFLAGS = 

if COND_GCOV
AM_CFLAGS = --coverage
AM_LDFLAGS += -lgcov
endif # COND_GCOV

GCOV_FILES = *.gcda *.gcno *.gcov

CLEANFILES = $(GCOV_FILES) *~

clean-local:
	-rm -rf *.dSYM


EXTRA_DIST = lexer.l parser.y

bin_PROGRAMS = freesasa
lib_LIBRARIES = libfreesasa.a
noinst_PROGRAMS = example
include_HEADERS = freesasa.h
libfreesasa_a_SOURCES = classifier.c classifier.h protor.c oons.c \
	coord.c coord.h pdb.c pdb.h \
	sasa_lr.c sasa_sr.c structure.c \
	freesasa.c freesasa.h freesasa_internal.h \
	nb.h nb.c util.c \
	selection.h selection.c $(lp_output)
freesasa_SOURCES = main.c 
example_SOURCES = example.c
freesasa_LDADD = libfreesasa.a
example_LDADD = libfreesasa.a

lp_output = lexer.c lexer.h parser.c parser.h

if GENERATE_PARSER
$(lp_output): lexer.l parser.y
	@LEX@ lexer.l
	@YACC@ parser.y

CLEANFILES += $(lp_output)

endif # GENERATE_PARSER


