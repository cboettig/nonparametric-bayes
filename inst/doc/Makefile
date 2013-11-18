SOURCEFILE = nonparametric-bayes



all: 
	make pdf


github:  ${SOURCEFILE}.Rmd ${SOURCEFILE}.bib
	./knit github
	pandoc -s -S --biblio ${SOURCEFILE}.bib ${SOURCEFILE}.md -o ${SOURCEFILE}_github.md
	pandoc -s --template=rstudio-headers.html ${SOURCEFILE}_github.md -o  ${SOURCEFILE}.html

md: ${SOURCEFILE}.Rmd ${SOURCEFILE}.bib
	pandoc -s -S --biblio ${SOURCEFILE}.bib ${SOURCEFILE}.md -o ${SOURCEFILE}.md

odt: ${SOURCEFILE}.Rmd ${SOURCEFILE}.bib
	./knit docx
	pandoc -s -S --biblio ${SOURCEFILE}.bib ${SOURCEFILE}.md --reference-odt=reference.odt -o ${SOURCEFILE}.odt


docx:  ${SOURCEFILE}.Rmd ${SOURCEFILE}.bib
	./knit docx
	pandoc -s -S --biblio ${SOURCEFILE}.bib ${SOURCEFILE}.md --reference-docx=reference.docx -o ${SOURCEFILE}.docx


pdf: ${SOURCEFILE}.Rmd ${SOURCEFILE}.bib
	./knit pdf
	pandoc -s -S --biblio ${SOURCEFILE}.bib ${SOURCEFILE}.md --template=elsarticle.latex -o ${SOURCEFILE}.tex
	xelatex ${SOURCEFILE}.tex
	xelatex ${SOURCEFILE}.tex
	xelatex ${SOURCEFILE}.tex

epub: ${SOURCEFILE}.Rmd ${SOURCEFILE}.bib
	./knit epub
	pandoc -s -S --biblio ${SOURCEFILE}.bib ${SOURCEFILE}.md -o ${SOURCEFILE}.epub

bib: 
	grep -o "@\w*" ${SOURCEFILE}.Rmd | tr '\n' ',' | tr -d "@"  > refs.tmp  # Extract cited references
	cp auxtemplate refs.aux          ## Create the aux file (in 3 steps)
	cat refs.tmp >> refs.aux # need to avoid newline
	echo "}" >> refs.aux	
	bibtool -s -d --preserve.key.case=on -x refs.aux ~/Documents/Mendeley/bib/*.bib > ${SOURCEFILE}.bib

rebib:
	bibtool -s -d --preserve.key.case=on -x refs.aux ~/Documents/Mendeley/bib/*.bib > ${SOURCEFILE}.bib

clean:
	rm -rf *.log *.bbl *.blg *.out *.spl
