.PHONY: slides
pdf:
	pdflatex main.tex

slides:
	mkdir slides -p
	pdftocairo -png -r 300 main.pdf slides/main
