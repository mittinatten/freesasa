fig_sources = lnr_circles.tex lnr_slice.tex
gen_eps = ${subst .tex,.eps,${fig_sources}}
other_eps = time.eps precision.eps
fig_eps = $(gen_eps) $(other_eps)
fig_pdf = ${subst .eps,.pdf,${fig_eps}}

if HAVE_LATEX
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
pdf-local: $(fig_pdf) 
%.pdf : %.eps
	epstopdf $< -o $@ 
else
nopdf:
	@echo "epstopdf needed to generate pdf-figures"
pdf-local: nopdf
endif

if HAVE_LATEX 
if HAVE_PDF
  dist_doc_DATA = $(fig_eps) $(fig_pdf)
else 
  dist_doc_DATA = $(fig_eps)
endif 
else !HAVE_LATEX
dist_doc_DATA = $(other_eps)
endif !HAVE_LATEX

clean-local:
	rm -f *.dvi *.aux *.log *~ *.pdf
	rm $(gen_eps)