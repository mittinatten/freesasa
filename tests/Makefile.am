if USE_CHECK

if COND_GCOV
AM_CFLAGS = --coverage
AM_LDFLAGS = -lgcov
GCOV_FILES = *.gcda *.gcno *.gcov
else
GCOV_FILES = 
endif

EXTRA_DIST = data


TESTS = test
check_PROGRAMS = test
test_SOURCES = test.c test_pdb.c test_freesasa.c test_structure.c \
	test_classify.c test_coord.c test_adjacency.c

test_CFLAGS = -I$(top_srcdir)/src $(AM_CFLAGS) \
	-DDATADIR=\"$(top_srcdir)/tests/data/\"

test_LDADD = ../src/libfreesasa.a @CHECK_LIBS@

CLEANFILES = tmp/*  $(GCOV_FILES)

endif
