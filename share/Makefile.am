AM_LDFLAGS = 

if COND_GCOV
AM_CFLAGS = --coverage
AM_LDFLAGS += -lgcov
endif # COND_GCOV

GCOV_FILES = *.gcda *.gcno *.gcov

CLEANFILES = $(GCOV_FILES) *~

clean-local:
	-rm -rf *.dSYM


dataroot_DATA = dssp.config naccess.config oons.config protor.config


