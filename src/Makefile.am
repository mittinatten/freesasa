AM_LDFLAGS =
AM_CFLAGS =
freesasa_LDADD = -lstdc++
AM_CPPFLAGS = -I$(top_srcdir)/third-party/gemmi/include/

if COND_GCOV
AM_CFLAGS += --coverage
endif # COND_GCOV

GCOV_FILES = *.gcda *.gcno *.gcov

CLEANFILES = $(GCOV_FILES) *~ .deps/*

clean-local:
	-rm -rf *.dSYM

EXTRA_DIST = lexer.l parser.y

bin_PROGRAMS = freesasa
lib_LIBRARIES = libfreesasa.a
noinst_PROGRAMS = example
include_HEADERS = freesasa.h
libfreesasa_a_SOURCES = classifier.c classifier.h \
	classifier_protor.c classifier_oons.c classifier_naccess.c \
	coord.c coord.h pdb.c pdb.h log.c \
	sasa_lr.c sasa_sr.c structure.c node.c \
	freesasa.c freesasa.h freesasa_internal.h \
	nb.h nb.c util.c rsa.c \
	selection.h selection.c $(lp_output)
freesasa_SOURCES = main.cc cif.cc cif.hh
example_SOURCES = example.c
freesasa_LDADD += libfreesasa.a
example_LDADD = libfreesasa.a

lp_output = lexer.c lexer.h parser.c parser.h

if USE_JSON
libfreesasa_a_SOURCES += json.c
freesasa_LDADD += -ljson-c
example_LDADD += -ljson-c
endif # USE_JSON

if USE_XML
libfreesasa_a_SOURCES += xml.c
AM_CFLAGS += ${libxml2_CFLAGS}
freesasa_LDADD += ${libxml2_LIBS}
example_LDADD += ${libxml2_LIBS}
endif # USE_XML

if GENERATE_PARSER
$(lp_output): lexer.l parser.y
	@LEX@ --nounistd lexer.l
	@YACC@ parser.y

CLEANFILES += $(lp_output)

endif # GENERATE_PARSER
