IN=doc.tex
OUT=BFGDocumentation
DEL = .log .out .aux .toc .bbl .blg .run.xml -blx.bib

cleanfunc = for i in $(DEL); do rm -f $(OUT)$$i; done

all: bibliography

clean:
	$(cleanfunc)

simple: 
	$(cleanfunc)
	pdflatex --jobname=$(OUT) $(IN) 
	$(cleanfunc)


bibliography:
	$(cleanfunc)
	pdflatex --jobname=$(OUT) $(IN)
	bibtex  $(OUT).aux
	pdflatex --jobname=$(OUT) $(IN)
	$(cleanfunc)
