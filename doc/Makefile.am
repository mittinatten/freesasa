SUBDIRS = fig
manual = manual.tex
fig_eps = fig/lnr_circles.eps fig/lnr_slice.eps fig/precision.eps fig/threads.eps
fig_pdf = $(subst .eps,.pdf,$(fig_eps))

CLEANFILES =
ALLFILES =

EXTRA_DIST = $(manual) doxy_main.md GPL.md README.md doxygen-footer.html doxygen-header.html

if HAVE_LATEX

doc_DATA = manual.pdf manual.ps $(fig_eps) $(fig_pdf)

$(fig_pdf):
	@cd fig
	@make pdf
	@cd ..

$(fig_eps): 
	@cd fig
	@make ps
	@cd ..

manual.dvi: $(srcdir)/$(manual) $(fig_eps)
	latex -interaction=batchmode $<
	latex -interaction=batchmode $<

manual.ps: manual.dvi
	dvips -Ppdf manual.dvi -q

manual.pdf: $(srcdir)/$(manual)  $(fig_eps)
	pdflatex -interaction=batchmode $<
	pdflatex -interaction=batchmode $<
	@if [ ! -d "html/" ]; then mkdir html; fi
	cp manual.pdf html/

dvi-local: manual.dvi
ps-local: manual.ps
pdf-local: manual.pdf
ALLFILES += manual.dvi manual.ps manual.pdf

if BUILD_HTML_MANUAL

manual.xml: $(srcdir)/$(manual) $(fig_eps)
	latexml --quiet $< > manual.xml

html/manual.html: manual.xml
	cd html; latexmlpost ../manual.xml  --format=html5 --destination=tmp.html --nodefaultresources --navigationtoc=context 
	perl prepareHTML.pl html/tmp.html > $@
	rm html/tmp.html html/LaTeXML.*

html-local: html/manual.html

ALLFILES += html/manual.html
CLEANFILES += manual.xml

else # BUILD_HTML_MANUAL

html-local:
	@echo latexml not available, html manual will not be built.

endif # BUILD_HTML_MANUAL

else # HAVE_LATEX
nolatex:
	@echo LaTeX disabled, run configure script with option --enable-latex to enable
dvi-local: nolatex
ps-local: nolatex
pdf-local: nolatex


endif # HAVE_LATEX

if HAVE_DOXYGEN

doxyfile.stamp: $(top_srcdir)/src/freesasa.h doxy_main.md $(top_srcdir)/README.md Doxyfile \
	doxygen-footer.html doxygen-header.html
	@if [ ! -d "html/" ]; then mkdir html; fi
	$(DOXYGEN) Doxyfile
	echo Timestamp > doxyfile.stamp

CLEANFILES += doxyfile.stamp
ALLFILES += doxyfile.stamp

doxygen: doxyfile.stamp
else
doxygen:
	@echo Doxygen disabled, run configure-script with option --enable-doxygen to enable
endif

if BUILD_DOC
all-local: $(ALLFILES)
else
all-local:
endif

clean-local:
	rm -f *.ps *.dvi *.aux *.log *.pdf *.toc *~ 
	rm -rf html/doxygen
	rm -rf html/fig

