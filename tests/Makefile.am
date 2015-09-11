
EXTRA_DIST = data test-cli.sh.in run-test-cli.sh
check_PROGRAMS =
TESTS = 

if COND_GCOV
AM_CFLAGS = --coverage
AM_LDFLAGS = -lgcov
endif

GCOV_FILES = *.gcda *.gcno *.gcov

if USE_CHECK

TESTS += libtest statictest
check_PROGRAMS += libtest statictest
libtest_SOURCES = libtest.c test_pdb.c test_freesasa.c test_structure.c \
	test_classify.c test_coord.c test_verlet.c
statictest_SOURCES = statictest.c test_verlet_static.c

libtest_CFLAGS = -I$(top_srcdir)/src $(AM_CFLAGS) \
	-DDATADIR=\"$(top_srcdir)/tests/data/\"
statictest_CFLAGS = $(libtest_CFLAGS)

libtest_LDADD = ../src/libfreesasa.a -lcheck
statictest_LDADD = -lcheck

endif 

if RUN_CLI_TESTS # on by default
TESTS += run-test-cli.sh
endif

CLEANFILES = tmp/*  $(GCOV_FILES) *~

