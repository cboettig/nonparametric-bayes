SOURCEFILE = nonparametric-bayes

make:
#	make refs.aux 
	make ${SOURCEFILE}.bib
	make ${SOURCEFILE}.tex
	make ${SOURCEFILE}.pdf


${SOURCEFILE}.tex:  ${SOURCEFILE}.Rmd ${SOURCEFILE}.bib elsarticle.latex 
	./knit pdf
	pandoc -s -S --biblio ${SOURCEFILE}.bib ${SOURCEFILE}.md --template=elsarticle.latex -o ${SOURCEFILE}.tex

${SOURCEFILE}.pdf: ${SOURCEFILE}.tex 
	xelatex ${SOURCEFILE}.tex
	xelatex ${SOURCEFILE}.tex
	xelatex ${SOURCEFILE}.tex


${SOURCEFILE}.html:  ${SOURCEFILE}.Rmd ${SOURCEFILE}.bib
	./knit github
	pandoc -s -S --biblio ${SOURCEFILE}.bib ${SOURCEFILE}.md -o ${SOURCEFILE}_github.md
	pandoc -s --template=rstudio-headers.html ${SOURCEFILE}_github.md -o  ${SOURCEFILE}.html


# Build a list of references appearing in the source markdown
refs.aux: ${SOURCEFILE}.Rmd 
	grep -o "@\w*" ${SOURCEFILE}.Rmd | tr '\n' ',' | tr -d "@"  > refs.tmp  # Extract cited references
	sed '1i\\\citation{' refs.tmp > refs.aux
	echo "}" >> refs.aux	

${SOURCEFILE}.bib: refs.aux
	bibtool -s -d --preserve.key.case=on -x refs.aux ~/Documents/Mendeley/bib/*.bib > ${SOURCEFILE}.bib

clean:
	rm -rf *.log *.bbl *.blg *.out *.spl



## Additional formats 

${SOURCEFILE}.odt: ${SOURCEFILE}.Rmd ${SOURCEFILE}.bib
	./knit docx
	pandoc -s -S --biblio ${SOURCEFILE}.bib ${SOURCEFILE}.md --reference-odt=reference.odt -o ${SOURCEFILE}.odt
${SOURCEFILE}.docx:  ${SOURCEFILE}.Rmd ${SOURCEFILE}.bib
	./knit docx
	pandoc -s -S --biblio ${SOURCEFILE}.bib ${SOURCEFILE}.md --reference-docx=reference.docx -o ${SOURCEFILE}.docx
${SOURCEFILE}.epub: ${SOURCEFILE}.Rmd ${SOURCEFILE}.bib
	./knit epub
	pandoc -s -S --biblio ${SOURCEFILE}.bib ${SOURCEFILE}.md -o ${SOURCEFILE}.epub

