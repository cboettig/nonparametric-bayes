FROM rocker/hadleyverse
MAINTAINER Carl Boettiger cboettig@ropensci.org
RUN apt-get update \
  && apt-get install -y r-cran-R2jags r-cran-cairo \
  && installGithub.r --deps TRUE \
    cboettig/cboettigR \
    cboettig/pdg_control \
    cboettig/nonparametric-bayes
RUN git clone https://github.com/cboettig/nonparametric-bayes.git \
  && chown -R docker:docker /nonparametric-bayes

## Set up a copy in the rstudio home directory so it is easy to find
RUN mkdir /home/rstudio \ 
&& cp -r /nonparametric-bayes/manuscripts /home/rstudio/ \
&& chown -R docker:docker /home/rstudio



# RUN Rscript -e 'rmarkdown::render("/nonparametric-bayes/manuscripts/manuscript.Rmd")'
# RUN Rscript -e 'rmarkdown::render("/nonparametric-bayes/manuscripts/supplement.Rmd")'

