FROM rocker/hadleyverse
MAINTAINER Carl Boettiger cboettig@ropensci.org
RUN apt-get update \
  && apt-get install -y r-cran-R2jags r-cran-cairo \
  && installGithub.r --deps TRUE \
    cboettig/cboettigR \
    cboettig/pdg_control \
    cboettig/nonparametric-bayes
RUN git clone https://github.com/cboettig/nonparametric-bayes.git /manuscript \
  && chown -R docker:docker /manuscript
# RUN Rscript -e 'rmarkdown::render("/manuscript/manuscript.Rmd")'
# RUN Rscript -e 'rmarkdown::render("/manuscript/supplement.Rmd")'

