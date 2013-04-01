SOURCEFILE = nonparametric-bayes
SOURCEFILE_B = beyond-mechanistic-models



all: 
	make bib
	make pdf


github:
	./knit github
	pandoc -s -S --biblio ${SOURCEFILE}.bib ${SOURCEFILE}.md -o ${SOURCEFILE}_github.md
	pandoc -s --template=rstudio-headers.html ${SOURCEFILE}_github.md -o  ${SOURCEFILE}.html

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
	grep -o "@\w*" ${SOURCEFILE}.Rmd | tr '\n' ',' | tr -d "@"  > refs.tmp  # Extract cited references
	cp auxtemplate refs.aux          ## Create the aux file (in 3 steps)
	cat refs.tmp >> refs.aux
	echo "}" >> refs.aux	
	bibtool -s -d --preserve.key.case=on -x refs.aux ~/Documents/Mendeley/bib/*.bib > ${SOURCEFILE}.bib

fileb:
	./knitb fileb
	pandoc -s -S --biblio ${SOURCEFILE}.bib ${SOURCEFILE_B}.md --template=elsarticle.latex -o ${SOURCEFILE_B}.tex
	pdflatex ${SOURCEFILE_B}.tex
	pdflatex ${SOURCEFILE_B}.tex
	pdflatex ${SOURCEFILE_B}.tex

clean:
	rm -rf *.log *.bbl *.blg *.out *.spl
