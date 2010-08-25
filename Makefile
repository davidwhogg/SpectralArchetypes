.SUFFIXES: .tex .dvi .ps .pdf

all: quasar-redshift.pdf binary-quasar.pdf lrg-template.pdf

%.pdf: %.tex
	pdflatex $<
	- bash -c " ( grep Rerun $*.log && pdflatex $< ) || echo noRerun "
	- bash -c " ( grep Rerun $*.log && pdflatex $< ) || echo noRerun "
