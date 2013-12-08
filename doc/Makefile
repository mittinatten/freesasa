all: ms fig

.PHONY: fig

ms: manual.pdf

manual.pdf: manual.tex
	latex manual.tex
	dvips -Ppdf manual.dvi
	ps2pdf manual.ps

fig: 
	$(MAKE) -C fig/