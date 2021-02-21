EXTRA_DIST = data test-cli.in
check_PROGRAMS =
TESTS =
AM_CFLAGS =
AM_LDFLAGS =

if COND_GCOV
AM_CFLAGS += --coverage
endif # COND_GCOV

GCOV_FILES = *.gcda *.gcno *.gcov

if USE_CHECK

TESTS += test-api
check_PROGRAMS += test-api
test_api_SOURCES = test_main.c test_pdb.c test_freesasa.c test_structure.c \
	test_classifier.c test_coord.c test_nb.c test_selection.c tools.h tools.c \
	test_node.c

AM_CFLAGS += -I$(top_srcdir)/src -DDATADIR=\"$(top_srcdir)/tests/data/\" -DSHAREDIR=\"$(top_srcdir)/share/\"

test_api_LDADD = ../src/libfreesasa.a

if USE_JSON
test_api_SOURCES += test_json.c
test_api_LDADD += -ljson-c
endif # USE_JSON

if USE_XML
test_api_SOURCES += test_xml.c
test_api_LDADD += ${libxml2_LIBS}
AM_CFLAGS += ${libxml2_CFLAGS}
endif # USE_XML

endif # USE_CHECK

if RUN_CLI_TESTS # on by default
TESTS += test-cli
endif # RUN_CLI_TESTS


CLEANFILES = tmp/*  $(GCOV_FILES) *~ .deps/*

clean-local:
	-rm -rf *.dSYM
