#FILENAME=alignCompress
#TEXFILES=alignCompress.tex
#PLOTFILES=p1.eps p2.eps p3.eps p4.eps p5.eps p6.eps p7.eps   r1.eps r2.eps r3.eps r4.eps

FILENAME=bioinformatics2003
TEXFILES=bioinformatics2003.tex
PLOTFILES=roc_uni_0.eps roc_blend.eps
#PLOTFILES=roc_uni.eps roc_0.eps roc_uni_0.eps roc_blend.eps
POP_CONV_SCRIPT=../testRuns/populationPlot.pl
POP_CONV_SCRIPT2=../testRuns/populationPlot2.pl
REL_CONV_SCRIPT=../testRuns/relatedPlot.pl

$(FILENAME).ps: $(TEXFILES) $(PLOTFILES) $(FILENAME).bbl
#$(FILENAME).ps: $(TEXFILES) $(PLOTFILES)
	@echo @@@@@@@@@@@@@@@@@@@@@@@@ Pass 1 of latex @@@@@@@@@@@@@@@@@@@@@@@
	latex $(FILENAME)
	@echo @@@@@@@@@@@@@@@@@@@@@@@@ Pass 2 of latex @@@@@@@@@@@@@@@@@@@@@@@
	latex $(FILENAME)
	@echo @@@@@@@@@@@@@@@@@@@@@@@@ dvips           @@@@@@@@@@@@@@@@@@@@@@@
	dvips $(FILENAME)

$(FILENAME).bbl: $(TEXFILES)
	latex $(FILENAME)
	bibtex $(FILENAME)

roc_uni.eps: ../testRuns/res.population.20020924_2.txt $(POP_CONV_SCRIPT)
	$(POP_CONV_SCRIPT) $< $@

roc_0.eps: ../testRuns/res.population.20020924.txt $(POP_CONV_SCRIPT)
	$(POP_CONV_SCRIPT) $< $@

roc_uni_0.eps: ../testRuns/res.population.20020925.txt $(POP_CONV_SCRIPT)
	$(POP_CONV_SCRIPT) $< $@

roc_blend.eps: ../testRuns/res.population2.20020920.txt $(POP_CONV_SCRIPT2)
	$(POP_CONV_SCRIPT2) $< $@




p1.eps: ../testRuns/res.population.20020913.txt $(POP_CONV_SCRIPT)
	$(POP_CONV_SCRIPT) $< $@

p2.eps: ../testRuns/res.population.20020915.txt $(POP_CONV_SCRIPT)
	$(POP_CONV_SCRIPT) $< $@

p3.eps: ../testRuns/res.population.20020913_2.txt $(POP_CONV_SCRIPT)
	$(POP_CONV_SCRIPT) $< $@

p4.eps: ../testRuns/res.population.20020914_2.txt $(POP_CONV_SCRIPT)
	$(POP_CONV_SCRIPT) $< $@

p5.eps: ../testRuns/res.population.20020915_2.txt $(POP_CONV_SCRIPT)
	$(POP_CONV_SCRIPT) $< $@

p6.eps: ../testRuns/res.population2.20020917.txt $(POP_CONV_SCRIPT2)
	$(POP_CONV_SCRIPT2) $< $@

p7.eps: ../testRuns/res.population2.20020920.txt $(POP_CONV_SCRIPT2)
	$(POP_CONV_SCRIPT2) $< $@


r1.eps: ../testRuns/res.related_20020827.txt $(REL_CONV_SCRIPT)
	$(REL_CONV_SCRIPT) $< $@

r2.eps: ../testRuns/res.related_20020827_2.txt $(REL_CONV_SCRIPT)
	$(REL_CONV_SCRIPT) $< $@

r3.eps: ../testRuns/res.related_20020902.txt $(REL_CONV_SCRIPT)
	$(REL_CONV_SCRIPT) $< $@

r4.eps: ../testRuns/res.related_20020903.txt $(REL_CONV_SCRIPT)
	$(REL_CONV_SCRIPT) $< $@

clean:
	rm -f $(FILENAME).ps *~ core *.aux *.log *.dvi *.toc *.bak *.lot *.lof *.bbl *.blg *.lol

realclean:
	rm -f $(FILENAME).ps $(PLOTFILES) *~ core *.aux *.log *.dvi *.toc *.bak *.lot *.lof *.bbl *.blg *.lol p?.eps r?.eps

