all:
	make manuscript.pdf
	make manuscript.tex

manuscript.md: manuscript.Rmd
	Rscript components/knit.R manuscript.Rmd	

manuscript.pdf: manuscript.md
	Rscript components/pandoc.R manuscript.md	

manuscript.tex: manuscript.md
	pandoc manuscript.md --template=components/elsarticle.latex --bibliography=components/references.bib --csl=components/ecology-letters.csl -o manuscript.tex

restore-cache:
	Rscript components/restore-cache.R

clear-cache:
	rm -rf components/cache
clean-reset:
	rm -rf manuscript.md manuscript.pdf components/cache components/figure
clean:


