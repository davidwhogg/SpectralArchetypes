LATEX       = pdflatex
BASH        = bash -c
ECHO        = echo
RM          = rm -rfv
RM_TMP      = ${RM} *.aux *.log *.toc
CHECK_RERUN = grep Rerun $*.log

all: archetypes.pdf

%.pdf: %.tex
	${LATEX} $<
	( ${CHECK_RERUN} && ${LATEX} $< ) || echo "Done."
	( ${CHECK_RERUN} && ${LATEX} $< ) || echo "Done."

clean:
	${RM_TMP}
