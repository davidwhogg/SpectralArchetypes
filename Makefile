.SUFFIXES: .tex .dvi .ps .pdf

all: archetypes.pdf

%.pdf: %.tex
	pdflatex $<
	pdflatex $<
	pdflatex $<
