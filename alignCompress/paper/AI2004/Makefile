#FILENAME=alignCompress
#TEXFILES=alignCompress.tex
#PLOTFILES=p1.eps p2.eps p3.eps p4.eps p5.eps p6.eps p7.eps   r1.eps r2.eps r3.eps r4.eps

FILENAME=AI2004
TEXFILES=AI2004.tex
EPS_PLOTFILES=roc_uni_0.eps roc_blend.eps
PDF_PLOTFILES=roc_uni_0.pdf roc_blend.pdf
#PLOTFILES=roc_uni.eps roc_0.eps roc_uni_0.eps roc_blend.eps

main: ps

pdf: $(TEXFILES) $(PDF_PLOTFILES) $(FILENAME).bbl
	@echo @@@@@@@@@@@@@@@@@@@@@@@@ Pass 1 of latex @@@@@@@@@@@@@@@@@@@@@@@
	pdflatex $(FILENAME)
	@echo @@@@@@@@@@@@@@@@@@@@@@@@ Pass 2 of latex @@@@@@@@@@@@@@@@@@@@@@@
	pdflatex $(FILENAME)
	@echo @@@@@@@@@@@@@@@@@@@@@@@@ Pass 3 of latex @@@@@@@@@@@@@@@@@@@@@@@
	pdflatex $(FILENAME)

ps: $(TEXFILES) $(EPS_PLOTFILES) $(FILENAME).bbl
	@echo @@@@@@@@@@@@@@@@@@@@@@@@ Pass 1 of latex @@@@@@@@@@@@@@@@@@@@@@@
	latex $(FILENAME)
	@echo @@@@@@@@@@@@@@@@@@@@@@@@ Pass 2 of latex @@@@@@@@@@@@@@@@@@@@@@@
	latex $(FILENAME)
	@echo @@@@@@@@@@@@@@@@@@@@@@@@ dvips           @@@@@@@@@@@@@@@@@@@@@@@
	dvips $(FILENAME)

$(FILENAME).bbl: $(TEXFILES)
	@echo @@@@@@@@@@@@@@@@@@@@@@@@ latex before bibtex @@@@@@@@@@@@@@@@@@@@@@@
	latex $(FILENAME)
	@echo @@@@@@@@@@@@@@@@@@@@@@@@ bibtex @@@@@@@@@@@@@@@@@@@@@@@
	bibtex $(FILENAME)

$(PDF_PLOTFILES): $(EPS_PLOTFILES)
	for F in $(EPS_PLOTFILES); do epstopdf $$F; done

clean:
	rm -f $(FILENAME).bbl $(FILENAME).dvi $(FILENAME).blg $(FILENAME).log $(FILENAME).aux
