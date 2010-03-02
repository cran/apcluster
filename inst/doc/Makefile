all: APCluster-Manual.pdf clean

clean:
	rm -f *.log *.aux *.bak *~ *.eps *.out *.toc *.blg *.tex Rplots.pdf
	rm -rf .Rd2dvi* APCluster-Manual-*.pdf

APCluster-Manual.pdf: APCluster-Manual.tex Makefile bioinf.cls bioinf-article.txi
	pdflatex APCluster-Manual
	pdflatex APCluster-Manual

APCluster-Manual.tex: APCluster-Manual.Rnw Makefile
	R CMD Sweave APCluster-Manual.Rnw
