all:
	make manuscript.pdf
	make supplement.pdf
	make ../README.md

install: components/install.R
	Rscript components/install.R

restore-cache:  
	make components/cache

manuscript.pdf: manuscript.Rmd
	Rscript -e "rmarkdown::render('manuscript.Rmd')"

supplement.pdf: supplement.Rmd
	Rscript -e "rmarkdown::render('manuscript.Rmd'); rmarkdown::render('supplement.Rmd')"

clear-cache:
	rm -rf components/cache

clean-reset:
	rm -rf manuscript.md manuscript.pdf components/cache components/figure

clean:
	rm -rf manuscript.md *.log *.spl *.aux *.out

# Non portable, just for me to update the source cache
deploy-cache:
	cd components; tar -pczf nonparametric-bayes-cache.tar.gz cache; scp nonparametric-bayes-cache.tar.gz ten:~/public_html/data/


# Only runs if cache has been cleared
components/cache: components/restore-cache.R
	Rscript components/restore-cache.R

