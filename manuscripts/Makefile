all:
	make manuscript.pdf

manuscript.md: manuscript.Rmd
	Rscript components/knit.R manuscript.Rmd	

manuscript.pdf: manuscript.md
	Rscript components/pandoc.R manuscript.md	

restore-cache:
	Rscript components/restore-cache.R

clear-cache:
	rm -rf components/cache
clean-reset:
	rm -rf manuscript.md manuscript.pdf components/cache components/figure
clean:


