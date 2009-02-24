.SUFFIXES: .tex .dvi .ps .pdf

all: archetypes.pdf emline.pdf

%.pdf: %.tex
	pdflatex $<
	pdflatex $<
	pdflatex $<
