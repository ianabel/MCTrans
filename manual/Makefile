
MCTrans.pdf: MCTrans.tex macros.tex references.bib gyrokinetics-macros.tex
	pdflatex MCTrans
	pdflatex MCTrans
	bibtex MCTrans
	pdflatex MCTrans
	bibtex MCTrans
	pdflatex MCTrans

clean:
	rm -f MCTrans.pdf MCTrans.out MCTrans.log MCTrans.bbl MCTrans.blg MCTrans.aux macros.aux gyrokinetic-macros.aux
