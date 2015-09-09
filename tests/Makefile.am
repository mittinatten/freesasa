
EXTRA_DIST = data test-cli.sh.in run-test-cli.sh
check_PROGRAMS =
TESTS = 

if COND_GCOV
AM_CFLAGS = --coverage
AM_LDFLAGS = -lgcov
endif

GCOV_FILES = *.gcda *.gcno *.gcov

if USE_CHECK

TESTS += libtest
check_PROGRAMS += libtest
libtest_SOURCES = libtest.c test_pdb.c test_freesasa.c test_structure.c \
	test_classify.c test_coord.c test_verlet.c

libtest_CFLAGS = -I$(top_srcdir)/src $(AM_CFLAGS) \
	-DDATADIR=\"$(top_srcdir)/tests/data/\"

libtest_LDADD = ../src/libfreesasa.a -lcheck

endif 

# we want this to run second if Check is active, since it is much slower
TESTS += run-test-cli.sh

CLEANFILES = tmp/*  $(GCOV_FILES) *~

