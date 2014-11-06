FROM rocker/hadleyverse
MAINTAINER Carl Boettiger cboettig@ropensci.org
RUN apt-get update \
  && apt-get install -y r-cran-R2jags \
  && installGithub.r --deps TRUE \
    cboettig/cboettigR \
    cboettig/pdg_control \
    cboettig/nonparametric-bayes
COPY . /manuscript

# RUN Rscript -e 'rmarkdown::render("manuscript/manuscript.Rmd")'
# RUN Rscript -e 'rmarkdown::render("manuscript/supplement.Rmd")'

