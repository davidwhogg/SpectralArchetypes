.SUFFIXES: .tex .ps .pdf .png .jpg

ALLPLOTS = $(patsubst %.ps,%.jpg,$(wildcard paper_plots/*.ps))
all: hmf_method.pdf

hmf_method.pdf: hmf_method.tex $(ALLPLOTS)

%.jpg: %.ps
	convert -density 150 $< $*.jpg

%.pdf: %.tex
	pdflatex $<
	- bash -c " ( grep Rerun $*.log && pdflatex $< ) || echo noRerun "
	- bash -c " ( grep Rerun $*.log && pdflatex $< ) || echo noRerun "
