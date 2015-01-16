FROM rocker/hadleyverse
MAINTAINER Carl Boettiger cboettig@ropensci.org
RUN apt-get update \
  && apt-get install -y r-cran-R2jags r-cran-cairo \
  && installGithub.r --deps TRUE \
    cboettig/cboettigR \
    cboettig/pdg_control \
    cboettig/nonparametric-bayes
COPY . /home/rstudio/pdg-control
WORKDIR /home/rstudio/pdg-control


# RUN Rscript -e 'rmarkdown::render("/home/rstudio/pdg-control/manuscript.Rmd")'
# RUN Rscript -e 'rmarkdown::render("/home/rstudio/pdg-control/supplement.Rmd")'

