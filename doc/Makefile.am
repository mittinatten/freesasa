CLEANFILES =
ALLFILES =

EXTRA_DIST = doxy-main.md license.md README.md \
	doxygen-footer.html doxygen-header.html \
	fig man/freesasa.1

man1_MANS = man/freesasa.1

if HAVE_DOXYGEN

doxyfile.stamp: $(top_srcdir)/src/freesasa.h $(top_srcdir)/src/example.c \
	doxy-main.md $(top_srcdir)/README.md Doxyfile \
	doxygen-footer.html doxygen-header.html
	@if [ ! -d "html/" ]; then mkdir html; fi
	@if [ ! -d "html/fig/" ]; then mkdir html/fig; fi
	cp fig/* html/fig/
	$(DOXYGEN) Doxyfile
	echo Timestamp > doxyfile.stamp

CLEANFILES += doxyfile.stamp
ALLFILES += doxyfile.stamp

doxygen: doxyfile.stamp

all-local: doxyfile.stamp

else
doxygen:
	@echo Doxygen disabled, run configure-script with option --enable-doxygen to enable
endif

clean-local:
	rm -f *.log *~
	rm -rf html/doxygen
	rm -rf html/fig
