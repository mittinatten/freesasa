SUBDIRS = src doc tests bindings

EXTRA_DIST = README.md share/protor.config share/naccess.config share/dssp.config share/oons.config \
	scripts/chemcomp2config.pl scripts/config2c.pl

# we want to test all features for dist-check
DISTCHECK_CONFIGURE_FLAGS = --enable-check

CLEANFILES = *~ share/*~ scripts/*~
