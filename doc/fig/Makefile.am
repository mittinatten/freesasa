fig_sources = lnr_circles.tex lnr_slice.tex
gen_eps = ${subst .tex,.eps,${fig_sources}}
other_eps = scaling.eps precision.eps threads.eps
fig_eps = $(gen_eps) $(other_eps)
fig_pdf = ${subst .eps,.pdf,${fig_eps}}
fig_svg = ${subst .pdf,.svg,${fig_pdf}}

EXTRA_DIST = $(fig_sources) $(other_eps)

if HAVE_LATEX

dist_doc_DATA = $(other_eps)

doc_DATA = $(fig_eps) $(fig_pdf) 

dvi-local: $(gen_eps)

%.dvi : %.tex
	latex -interaction=batchmode $< ${subst .tex,.dvi,$<}
%.eps : %.dvi
	dvips -E ${subst .tex,.dvi,$<} -o $@ -q

%.pdf : %.eps
	epstopdf $< -o $@ 

if BUILD_HTML_MANUAL

doc_DATA += $(fig_svg)

%.svg : %.pdf
	pdftocairo -svg $< 

pdf-local: $(fig_pdf) $(fig_svg)

else # BUILD_HTML_MANUAL

pdf-local: $(fig_pdf)

endif # BUILD_HTML_MANUAL

else # HAVE_LATEX

nolatex:
#	@echo "LaTeX needed to generate figures for manual"
dvi-local: nolatex
pdf-local: nolatex

endif # HAVE_LATEX

ps-local: dvi-local
all-local: dvi-local pdf-local

clean-local:
	rm -f *.dvi *.aux *.log *~ *.pdf *.svg
	rm -f $(gen_eps)
