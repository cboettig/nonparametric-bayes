MANUSCRIPT=manuscript


all:
	make render
	make ../README.md

install: components/install.R
	Rscript components/install.R
restore-cache:  
	make components/cache

render: 
	Rscript -e "rmarkdown::render('$(MANUSCRIPT).Rmd')"

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

