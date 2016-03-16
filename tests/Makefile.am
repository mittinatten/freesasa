EXTRA_DIST = data test-cli.in whole_lib_one_file.c
check_PROGRAMS =
TESTS = 
AM_CFLAGS =
AM_LDFLAGS =

if COND_GCOV
AM_CFLAGS += --coverage
AM_LDFLAGS += -lgcov
endif # COND_GCOV

GCOV_FILES = *.gcda *.gcno *.gcov

if USE_CHECK

TESTS += test-api test-static test-memerr
check_PROGRAMS += test-api test-static test-memerr
test_api_SOURCES = libtest.c test_pdb.c test_freesasa.c test_structure.c \
	test_classifier.c test_coord.c test_nb.c test_selection.c tools.h tools.c
test_static_SOURCES = test_static.c tools.h tools.c
test_memerr_SOURCES = test_memerr.c tools.h tools.c

AM_CFLAGS += -I$(top_srcdir)/src -DDATADIR=\"$(top_srcdir)/tests/data/\" -DSHAREDIR=\"$(top_srcdir)/share/\"

test_api_LDADD = ../src/libfreesasa.a -lcheck 
test_static_LDADD = ../src/parser.o ../src/lexer.o -lcheck 
test_memerr_LDADD = ../src/parser.o ../src/lexer.o -lcheck 

endif # USE_CHECK

if RUN_CLI_TESTS # on by default
TESTS += test-cli
endif # RUN_CLI_TESTS

CLEANFILES = tmp/*  $(GCOV_FILES) *~

clean-local:
	-rm -rf *.dSYM
