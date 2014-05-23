all:
	make install
	make restore-cache
	make manuscript.pdf

install: 
	Rscript components/install.R

manuscript.md: manuscript.Rmd
	Rscript components/knit.R manuscript.Rmd	

# Also makes manuscript.tex
manuscript.pdf: manuscript.md
	Rscript components/pandoc.R manuscript.md	

restore-cache:
	make components/cache

clear-cache:
	rm -rf components/cache

clean-reset:
	rm -rf manuscript.md manuscript.pdf components/cache components/figure

clean:
	rm -rf manuscript.md

# Non portable, just for me to update the source cache
deploy-cache:
	cd components; tar -pczf nonparametric-bayes-cache.tar.gz cache; scp nonparametric-bayes-cache.tar.gz ten:~/public_html/data/




# Only runs if cache has been cleared
components/cache: components/restore-cache.R
	Rscript components/restore-cache.R

