if COND_GCOV
AM_CFLAGS = --coverage
AM_LDFLAGS = -lgcov
GCOV_FILES = *.gcda *.gcno *.gcov
else
GCOV_FILES = 
endif

TESTS = test
check_PROGRAMS = test
test_SOURCES = test.c test_pdb.c test_sasalib.c test_structure.c \
	test_classify.c test_coord.c

test_CFLAGS = -I$(top_srcdir)/src $(AM_CFLAGS)

test_LDADD = ../src/libsasa.a @CHECK_LIBS@

CLEANFILES = tmp/*  $(GCOV_FILES)