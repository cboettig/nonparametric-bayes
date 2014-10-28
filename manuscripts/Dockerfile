FROM rocker/hadleyverse

RUN apt-get update \
  && apt-get install -y r-cran-R2jags \
  && install2.r --error \
    coda \
    emdbook \
    kernlab \
    MASS \
    mcmc \
    MCMCpack \
    modeest \
    pander \
    R2jags

RUN installGithub.r \
    cboettig/pdg_control \
    cboettig/nonparametric-bayes \
    cboettig/cboettigR


## Don't do this via the dockerfile, we want to be able to link the directory

#WORKDIR /home/docker/manuscript/
#COPY manuscript.Rmd /home/docker/manuscript/manuscript.Rmd
#COPY components/ /home/docker/manuscript/components/
#RUN Rscript -e 'rmarkdown::render("manuscript.Rmd")'
#RUN Rscript -e 'rmarkdown::render("supplement.Rmd")'

