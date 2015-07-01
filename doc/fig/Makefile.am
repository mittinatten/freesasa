fig_sources = lnr_circles.tex lnr_slice.tex
gen_eps = ${subst .tex,.eps,${fig_sources}}
other_eps = time.eps precision.eps
fig_eps = $(gen_eps) $(other_eps)
fig_pdf = ${subst .eps,.pdf,${fig_eps}}

dist_doc_DATA = time.eps precision.eps lnr_circles.tex lnr_slice.tex

if HAVE_LATEX

dist_doc_DATA += lnr_circles.eps lnr_slice.eps

dvi-local: $(gen_eps)
%.eps : %.tex
	latex -interaction=batchmode $<
	dvips -E ${subst .tex,.dvi,$<} -o $@ -q
else
nolatex:
	@echo "LaTeX needed to generate figures for manual"
dvi-local: nolatex
endif

ps-local: dvi-local

if HAVE_PDF

dist_doc_DATA += lnr_circles.pdf lnr_slice.pdf

pdf-local: $(fig_pdf) 
%.pdf : %.eps
	epstopdf $< -o $@ 
else
nopdf:
	@echo "epstopdf needed to generate pdf-figures"
pdf-local: nopdf
endif

clean-local:
	rm -f *.dvi *.aux *.log *~ *.pdf
	rm -f $(gen_eps)
