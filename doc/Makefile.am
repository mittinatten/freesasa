docfiles = manual.tex

if HAVE_LATEX
doc_DATA = manual.pdf
MANNAME = manual
MANTEXSRC = $(MANNAME).tex
MANAUX = $(MANNAME).aux
MANPDF = $(MANNAME).pdf
MANPS = $(MANNAME).ps
MANDVI = $(MANNAME).dvi

CLEANFILES = $(MANPDF) $(MANPS) $(MANDVI) $(MANNAME).log $(MANNAME).idx \
	$(MANNAME).toc $(MANAUX)


$(MANPDF): $(MANTEXSRC)
	latex manual.tex
	dvips -Ppdf manual.dvi
	ps2pdf manual.ps
endif

dist_doc_DATA = $(docfiles)