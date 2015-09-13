
EXTRA_DIST = data test-cli.in run-test-cli whole_lib_one_file.c
check_PROGRAMS =
TESTS = 

if COND_GCOV
AM_CFLAGS = --coverage
AM_LDFLAGS = -lgcov
endif

GCOV_FILES = *.gcda *.gcno *.gcov

if USE_CHECK

static_tests = static-nb static-sasa-lr static-structure 

TESTS += libtest $(static_tests) broken-malloc
check_PROGRAMS += libtest $(static_tests) broken-malloc
libtest_SOURCES = libtest.c test_pdb.c test_freesasa.c test_structure.c \
	test_classify.c test_coord.c test_nb.c
static_nb_SOURCES = static_nb.c
static_sasa_lr_SOURCES = static_sasa_lr.c
static_structure_SOURCES = static_structure.c
broken_malloc_SOURCES = broken_malloc.c

CFLAGS = -I$(top_srcdir)/src $(AM_CFLAGS) \
	-DDATADIR=\"$(top_srcdir)/tests/data/\"

LDADD = -lcheck
libtest_LDADD = ../src/libfreesasa.a -lcheck

endif 

if RUN_CLI_TESTS # on by default
TESTS += run-test-cli
endif

CLEANFILES = tmp/*  $(GCOV_FILES) *~

