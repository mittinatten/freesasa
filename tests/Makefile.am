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

TESTS += libtest static-test broken-malloc
check_PROGRAMS += libtest static-test broken-malloc
libtest_SOURCES = libtest.c test_pdb.c test_freesasa.c test_structure.c \
	test_classifier.c test_coord.c test_nb.c test_selection.c tools.h tools.c
static_test_SOURCES = test_static.c tools.h
broken_malloc_SOURCES = broken_malloc.c tools.h 

AM_CFLAGS += -I$(top_srcdir)/src -DDATADIR=\"$(top_srcdir)/tests/data/\" -DSHAREDIR=\"$(top_srcdir)/share/\"

libtest_LDADD = ../src/libfreesasa.a -lcheck
static_test_LDADD = ../src/parser.o ../src/lexer.o -lcheck 
broken_malloc_LDADD = ../src/parser.o ../src/lexer.o -lcheck 

endif # USE_CHECK

if RUN_CLI_TESTS # on by default
TESTS += test-cli
endif # RUN_CLI_TESTS

CLEANFILES = tmp/*  $(GCOV_FILES) *~

