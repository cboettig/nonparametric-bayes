SOURCEFILE = nonparametric-bayes



all: 
	make bib
	make pdf


github: ${SOURCEFILE}.html
	./knit github
	pandoc -s -S --biblio ${SOURCEFILE}.bib ${SOURCEFILE}.md -o ${SOURCEFILE}_github.md
	pandoc -s --template=rstudio-headers.html ${SOURCEFILE}_github.md -o  ${SOURCEFILE}.html

md:
	pandoc -s -S --biblio ${SOURCEFILE}.bib ${SOURCEFILE}.md -o ${SOURCEFILE}.md

odt: ${SOURCEFILE}.odt
	./knit docx
	pandoc -s -S --biblio ${SOURCEFILE}.bib ${SOURCEFILE}.md --reference-odt=reference.odt -o ${SOURCEFILE}.odt


docx:  ${SOURCEFILE}.docx
	./knit docx
	pandoc -s -S --biblio ${SOURCEFILE}.bib ${SOURCEFILE}.md --reference-docx=reference.docx -o ${SOURCEFILE}.docx


pdf: ${SOURCEFILE}.pdf
	./knit pdf
	pandoc -s -S --biblio ${SOURCEFILE}.bib ${SOURCEFILE}.md --template=elsarticle.latex -o ${SOURCEFILE}.tex
	pdflatex ${SOURCEFILE}.tex
	pdflatex ${SOURCEFILE}.tex
	pdflatex ${SOURCEFILE}.tex

epub: ${SOURCEFILE}.epub
	./knit epub
	pandoc -s -S --biblio ${SOURCEFILE}.bib ${SOURCEFILE}.md -o ${SOURCEFILE}.epub

bib: ${SOURCEFILE}.bib
	grep -o "@\w*" ${SOURCEFILE}.Rmd | tr '\n' ',' | tr -d "@"  > refs.tmp  # Extract cited references
	cp auxtemplate refs.aux          ## Create the aux file (in 3 steps)
	cat refs.tmp >> refs.aux
	echo "}" >> refs.aux	
	bibtool -s -d --preserve.key.case=on -x refs.aux ~/Documents/Mendeley/bib/*.bib > ${SOURCEFILE}.bib


clean:
	rm -rf *.log *.bbl *.blg *.out *.spl
