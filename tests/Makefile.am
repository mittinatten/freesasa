
EXTRA_DIST = data test-cli.sh.in run-test-cli.sh
check_PROGRAMS =
TESTS = 

if COND_GCOV
AM_CFLAGS = --coverage
AM_LDFLAGS = -lgcov
endif

GCOV_FILES = *.gcda *.gcno *.gcov

if USE_CHECK

static_tests = static-verlet static-sasa-lr

TESTS += libtest $(static_tests)
check_PROGRAMS += libtest $(static_tests)
libtest_SOURCES = libtest.c test_pdb.c test_freesasa.c test_structure.c \
	test_classify.c test_coord.c test_verlet.c
static_verlet_SOURCES = static_verlet.c
static_sasa_lr_SOURCES = static_sasa_lr.c

CFLAGS = -I$(top_srcdir)/src $(AM_CFLAGS) \
	-DDATADIR=\"$(top_srcdir)/tests/data/\"

LDADD = -lcheck
libtest_LDADD = ../src/libfreesasa.a -lcheck

endif 

if RUN_CLI_TESTS # on by default
TESTS += run-test-cli.sh
endif

CLEANFILES = tmp/*  $(GCOV_FILES) *~

