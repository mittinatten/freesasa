SUBDIRS = fig
manual = manual.tex
fig_eps = fig/lnr_circles.eps fig/lnr_slice.eps fig/time.eps fig/precision.eps
fig_pdf = $(subst .eps,.pdf,$(fig_eps))

EXTRA_DIST = $(manual)
dist_doc_DATA = 

if HAVE_LATEX

dist_doc_DATA += manual.dvi manual.ps

$(fig_eps): 
	@cd fig
	@make ps
	@cd ..

manual.ps: manual.dvi
	dvips -Ppdf manual.dvi -q
manual.dvi: $(manual) $(fig_eps)
	latex -interaction=batchmode manual.tex
	latex -interaction=batchmode manual.tex


dvi-local: manual.dvi

ps-local: manual.ps

else
nolatex:
	@echo "LaTeX needed to build .dvi and .ps manual"
dvi-local: nolatex
ps-local: nolatex
endif

if HAVE_PDF

#dist_doc_DATA += manual.pdf

$(fig_pdf):
	@cd fig
	@make pdf
	@cd ..

pdf-local: manual.pdf 
manual.pdf: $(manual) $(fig_pdf)
	pdflatex -interaction=batchmode manual.tex 
	pdflatex -interaction=batchmode manual.tex 
else
pdf-local: 
	@echo "pdflatex needed to build .pdf-manual"
endif


if HAVE_DOXYGEN

doxyfile.stamp: $(top_srcdir)/src/*.h
	$(DOXYGEN) Doxyfile
	echo Timestamp > doxyfile.stamp

CLEANFILES = doxyfile.stamp 

all-local: doxyfile.stamp
else
all-local:	
	@echo "Doxygen needed to build full reference manual."
endif

clean-local:
	rm -f *.ps *.dvi *.aux *.log *.pdf *.toc *~ 
	rm -rf doxygen/
