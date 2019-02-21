Sources for: Avoiding tipping points in fisheries management through Gaussian Process Dynamic Programming
========================================================================================================

[![Build Status](http://server.carlboettiger.info:88/api/badge/github.com/cboettig/nonparametric-bayes/status.svg?branch=master)](http://server.carlboettiger.info:88/github.com/cboettig/nonparametric-bayes)

<!--
Manuscript code is run automatically any time changes are made to
this repository.  The badge above summarizes the status and contains
links to further details. 
--> 

This directory contains the `.Rmd` files for the manuscript and
supplement. These files use functions from the R package from this
repository (`nonparametric-bayes`) but also several additional R packages
and other software which are responsible for the parametric comparisons,
as well as software responsible merely for the formatting to go from
`.Rmd` to a `.pdf` file (e.g. `knitr`, `pandoc`, `LaTeX`).

All of these additional dependencies means that rebuilding the
manuscript from scratch means more software to install and more
points of failure. To help you get started, we provide a copy
of the software environment we use in this paper as a Docker
container. If you'd rather not worry about the installation details,
[Install Docker](https://docs.docker.com/installation)
on your laptop or a [user-friendly cloud
server](https://www.digitalocean.com/?refcode=08c6ac401e49) [^1]

The badge above is automatically generated when the manuscript
is rebuilt on a remote server. This provides a convenient way
to confirm that everything is still working.


Docker Quickstart
-----------------

[Install Docker](https://docs.docker.com/installation) on your laptop or server. Then from the `manuscripts/` directory run: 

```bash
docker build -t cboettig/nonparametric-bayes .
```


```
docker run --rm -ti -v $(pwd):/data cboettig/nonparametric-bayes R -e "rmarkdown::render('manuscript.Rmd')"
```


You can now open `manuscript.Rmd` or `supplement.Rmd` files in RStudio and run the code interactively or compile the pdfs from scratch.  See [RStudio's rmarkdown](http://rmarkdown.rstudio.com/) for details. 


Local Installation requirements
--------------------------------

- RStudio (includes `pandoc`)
- A LaTeX environment
- The `nonparametric-bayes` package [from github](http://github.com/cboettig/nonparametric-bayes), including all packages in the [DESCRIPTION](../DESCRIPTION) file (by default suggested packages are not installed). 

--------

[1]: this link provides $10 credit to Digital Ocean; enough to run a 1GB instance for a month.  Disclosure: If you subsequently become a Digital Ocean customer the authors receive a credit for the referral. 


