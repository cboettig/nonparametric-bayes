SOURCEFILE = nonparametric-bayes



all: 

github:
	./knit github
	pandoc -s -S --biblio ${SOURCEFILE}.bib ${SOURCEFILE}.md -o ${SOURCEFILE}_github.md

md:
	pandoc -s -S --biblio ${SOURCEFILE}.bib ${SOURCEFILE}.md -o ${SOURCEFILE}.md

odt:
	./knit docx
	pandoc -s -S --biblio ${SOURCEFILE}.bib ${SOURCEFILE}.md --reference-odt=reference.odt -o ${SOURCEFILE}.odt


docx:
	./knit docx
	pandoc -s -S --biblio ${SOURCEFILE}.bib ${SOURCEFILE}.md --reference-docx=reference.docx -o ${SOURCEFILE}.docx


pdf:
	./knit pdf
	pandoc -s -S --biblio ${SOURCEFILE}.bib ${SOURCEFILE}.md --template=elsarticle.latex -o ${SOURCEFILE}.tex
	pdflatex ${SOURCEFILE}.tex
	pdflatex ${SOURCEFILE}.tex
	pdflatex ${SOURCEFILE}.tex

epub:
	./knit epub
	pandoc -s -S --biblio ${SOURCEFILE}.bib ${SOURCEFILE}.md -o ${SOURCEFILE}.epub

bib:
	bibtool -s -d --preserve.key.case=on -x citations.aux ~/Documents/Mendeley/bib/*.bib > ${SOURCEFILE}.bib


clean:
	rm -rf *.log *.bbl *.blg *.out *.spl
