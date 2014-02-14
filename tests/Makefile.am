TESTS = test
check_PROGRAMS = test
test_SOURCES = test.c test_pdb.c test_sasalib.c test_structure.c \
	test_classify.c
test_CFLAGS = -I$(top_srcdir)/src
test_LDADD = ../src/libsasa.a @CHECK_LIBS@

CLEANFILES = tmp/*