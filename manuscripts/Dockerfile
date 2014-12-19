FROM rocker/hadleyverse
MAINTAINER Carl Boettiger cboettig@ropensci.org
RUN apt-get update \
  && apt-get install -y r-cran-R2jags r-cran-cairo \
  && installGithub.r --deps TRUE \
    cboettig/cboettigR \
    cboettig/pdg_control \
    cboettig/nonparametric-bayes
COPY . /home/rstudio
RUN chown -R docker:docker /home/rstudio
WORKDIR /home/rstudio


# RUN Rscript -e 'rmarkdown::render("/home/rstudio/manuscript.Rmd")'
# RUN Rscript -e 'rmarkdown::render("/home/rstudio/supplement.Rmd")'

