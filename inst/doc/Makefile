LATEXFILES = *.aux\
        *.bbl\
        *.blg\
        *.ilg\
        *.log\
        *.nlo\
        *.nls\
        *.toc\
        *.aux\
	*.map\
	*.sh\
	sequences.tex\
        *.out

vignette:
	$(R_HOME)/bin/R CMD Sweave PWMEnrich.Rnw
	pdflatex PWMEnrich.tex
	bibtex PWMEnrich
	pdflatex PWMEnrich.tex
	pdflatex PWMEnrich.tex
	rm -f *-plots.pdf
	rm -f $(LATEXFILES)
	rm -f *~
clean:
	rm -f $(LATEXFILES)
	rm -f *~

