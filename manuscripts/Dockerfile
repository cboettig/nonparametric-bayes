FROM rocker/hadleyverse
MAINTAINER Carl Boettiger cboettig@ropensci.org

## Install software dependencies
RUN apt-get update \
  && apt-get install -y \
    r-cran-cairo \
    r-cran-R2jags \
  && installGithub.r --deps TRUE \
    cboettig/cboettigR \
    cboettig/pdg_control \
    cboettig/nonparametric-bayes

## Copy the manuscripts/ dir ('.' is relative to Dockerfile location) onto the container.
COPY . /home/rstudio/manuscripts
WORKDIR /home/rstudio/manuscripts

## Uncomment to build manuscript + cache as part of container build step.
# RUN Rscript -e 'rmarkdown::render("/home/rstudio/manuscripts/manuscript.Rmd")'
# RUN Rscript -e 'rmarkdown::render("/home/rstudio/manuscripts/supplement.Rmd")'
