
all: memory.pdf epilepsy.pdf
# all: attention.pdf phaseByTopic.pdf rhythms.pdf travelingWaves.pdf

epilepsy.pdf: epilepsy.dvi
epilepsy.dvi: epilepsy.tex epilepsy.bib

memory.pdf: memory.tex memory.bib

reactivations.pdf: reactivations.dvi
reactivations.dvi: reactivations.tex reactivations.bib

attention.pdf: attention.tex attention.bib

phaseByTopic.pdf: phaseByTopic.tex rhythms.bib pac.bib

rhythms.pdf: rhythms.tex rhythms.bib

travelingWaves.pdf: travelingWaves.tex travelingWaves.bib

%.dvi: %.tex
	latex $*.tex
	bibtex $*.aux
	latex $*.tex
	latex $*.tex

%.pdf: %.dvi
	dvipdf $*.dvi
