all: vblr.pdf

%.dvi: %.tex
	latex $<
	bibtex $*.aux
	latex $<
	latex $<

%.pdf: %.dvi
	dvipdf $<
