CLEANFILES =
ALLFILES =

EXTRA_DIST = doxy_main.md GPL.md README.md \
	doxygen-footer.html doxygen-header.html \
	HOWTO_BUILD_WEBPAGES fig/

if HAVE_DOXYGEN

doxyfile.stamp: $(top_srcdir)/src/freesasa.h doxy_main.md $(top_srcdir)/README.md Doxyfile \
	doxygen-footer.html doxygen-header.html $(top_srcdir)/bindings/python/freesasa.pyx
	@if [ ! -d "html/" ]; then mkdir html; fi
	@if [ ! -d "html/fig/" ]; then mkdir html/fig; fi
	cp fig/* html/fig/
	$(DOXYGEN) Doxyfile
	echo Timestamp > doxyfile.stamp

CLEANFILES += doxyfile.stamp
ALLFILES += doxyfile.stamp

doxygen: doxyfile.stamp
else
doxygen:
	@echo Doxygen disabled, run configure-script with option --enable-doxygen to enable
endif

all-local: doxyfile.stamp

clean-local:
	rm -f *.log *~ 
	rm -rf html/doxygen
	rm -rf html/fig

