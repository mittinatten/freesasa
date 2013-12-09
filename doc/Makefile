all: fig ms

.PHONY: fig ms

ms: manual.tex
	latex manual.tex
	dvips -Ppdf manual.dvi
	ps2pdf manual.ps

fig: 
	$(MAKE) -C fig/