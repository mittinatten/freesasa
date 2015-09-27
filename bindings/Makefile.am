TESTS =

cython_src = python/cfreesasa.pxd python/freesasa.pyx 

EXTRA_DIST = $(cython_src) python/setup.py python/test.py python/data test-python.sh

if HAVE_CYTHON

TESTS += test-python.sh

python_bindings = python/freesasa.so python/freesasa.c

pylib: $(cython_src) $(top_srcdir)/src/freesasa.h $(top_srcdir)/src/libfreesasa.a
	@cd $(top_srcdir)/bindings/python/; \
		@PYTHON@ setup.py build_ext --inplace 

install-exec-hook:
	@cd $(top_srcdir)/bindings/python/; \
		@PYTHON@ setup.py install

$(python_bindings): pylib

all-local: $(python_bindings)

endif #have-cython