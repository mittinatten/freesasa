SUBDIRS = fig
manual = manual.tex
fig_eps = fig/lnr_circles.eps fig/lnr_slice.eps fig/time.eps fig/precision.eps
fig_pdf = $(subst .eps,.pdf,$(fig_eps))

dist_doc_DATA = manual.pdf manual.ps manual.dvi

if HAVE_LATEX

manual.ps: manual.dvi
	dvips -Ppdf manual.dvi
manual.dvi: $(manual) $(fig_eps)
	latex manual.tex
	latex manual.tex

dvi-local: manual.dvi

ps-local: manual.ps

else
nolatex:
	@echo "LaTeX needed to build .dvi and .ps manual"
dvi-local: nolatex
ps-local: nolatex
endif

if HAVE_PDF
pdf-local: manual.pdf 
manual.pdf: $(manual) $(fig_pdf)
	pdflatex manual.tex
	pdflatex manual.tex
else
pdf-local: 
	@echo "pdflatec needed to build .pdf-manual"
endif

clean-local:
	rm -f *.ps *.dvi *.aux *.log *.pdf *.toc *~