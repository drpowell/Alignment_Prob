FILENAME=alignCompress
TEXFILES=alignCompress.tex
PLOTFILES=p1.eps p2.eps p3.eps p4.eps p5.eps r1.eps r2.eps r3.eps r4.eps
POP_CONV_SCRIPT=../testRuns/populationPlot.pl
REL_CONV_SCRIPT=../testRuns/relatedPlot.pl

#$(FILENAME).ps: $(TEXFILES) $(FILENAME).bbl
$(FILENAME).ps: $(TEXFILES) $(PLOTFILES)
	@echo @@@@@@@@@@@@@@@@@@@@@@@@ Pass 1 of latex @@@@@@@@@@@@@@@@@@@@@@@
	latex $(FILENAME)
	@echo @@@@@@@@@@@@@@@@@@@@@@@@ Pass 2 of latex @@@@@@@@@@@@@@@@@@@@@@@
	latex $(FILENAME)
	@echo @@@@@@@@@@@@@@@@@@@@@@@@ dvips           @@@@@@@@@@@@@@@@@@@@@@@
	dvips $(FILENAME)

$(FILENAME).bbl: $(TEXFILES)
	latex $(FILENAME)
	bibtex $(FILENAME)

p1.eps: ../testRuns/res.population.20020829_1.txt $(POP_CONV_SCRIPT)
	$(POP_CONV_SCRIPT) $< $@

p2.eps: ../testRuns/res.population.20020829_2.txt $(POP_CONV_SCRIPT)
	$(POP_CONV_SCRIPT) $< $@

p3.eps: ../testRuns/res.population.20020830.txt $(POP_CONV_SCRIPT)
	$(POP_CONV_SCRIPT) $< $@

p4.eps: ../testRuns/res.population.20020902.txt $(POP_CONV_SCRIPT)
	$(POP_CONV_SCRIPT) $< $@

p5.eps: ../testRuns/res.population.20020902_2.txt $(POP_CONV_SCRIPT)
	$(POP_CONV_SCRIPT) $< $@


r1.eps: ../testRuns/res.related_20020827.txt $(REL_CONV_SCRIPT)
	$(REL_CONV_SCRIPT) $< $@

r2.eps: ../testRuns/res.related_20020827_2.txt $(REL_CONV_SCRIPT)
	$(REL_CONV_SCRIPT) $< $@

r3.eps: ../testRuns/res.related_20020902.txt $(REL_CONV_SCRIPT)
	$(REL_CONV_SCRIPT) $< $@

r4.eps: ../testRuns/res.related_20020903.txt $(REL_CONV_SCRIPT)
	$(REL_CONV_SCRIPT) $< $@

clean:
	rm -f $(FILENAME).ps *~ core *.aux *.log *.dvi *.toc *.bak *.lot *.lof *.bbl *.blg *.lol p?.eps r?.eps

