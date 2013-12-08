src = $(wildcard *.tex)
eps = ${subst .tex,.eps,${src}}

all: ${eps}

%.eps : %.tex
	latex $<; dvips -E ${subst .tex,.dvi,$<} -o $@; 
