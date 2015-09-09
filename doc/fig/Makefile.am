fig_sources = lnr_circles.tex lnr_slice.tex
gen_eps = ${subst .tex,.eps,${fig_sources}}
other_eps = scaling.eps precision.eps threads.eps
fig_eps = $(gen_eps) $(other_eps)
fig_pdf = ${subst .eps,.pdf,${fig_eps}}

EXTRA_DIST = $(fig_sources) $(other_eps)

if HAVE_LATEX

dist_doc_DATA = $(other_eps)

doc_DATA = $(fig_eps) $(fig_pdf)

dvi-local: $(gen_eps)

%.dvi : %.tex
	latex -interaction=batchmode $< ${subst .tex,.dvi,$<}
%.eps : %.dvi
	dvips -E ${subst .tex,.dvi,$<} -o $@ -q

pdf-local: $(fig_pdf) 
%.pdf : %.eps
	epstopdf $< -o $@ 
else
nolatex:
#	@echo "LaTeX needed to generate figures for manual"
dvi-local: nolatex
pdf-local: nolatex
endif

ps-local: dvi-local
all-local: dvi-local pdf-local

clean-local:
	rm -f *.dvi *.aux *.log *~ *.pdf
	rm -f $(gen_eps)
