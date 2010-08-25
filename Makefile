.SUFFIXES: .tex .dvi .ps .pdf

all: quasar-redshift.pdf binary-quasar.pdf lrg-template.tex

%.pdf: %.tex
	pdflatex $<
	- bash -c " ( grep Rerun $*.log && pdflatex $< ) || echo noRerun "
	- bash -c " ( grep Rerun $*.log && pdflatex $< ) || echo noRerun "
