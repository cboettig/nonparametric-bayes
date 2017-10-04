FROM rocker/r-ver:3.2.5
MAINTAINER Carl Boettiger cboettig@ropensci.org

## Install software dependencies
RUN apt-get update \
  && apt-get install -y \
    libssl-dev \
    libxml2-dev \
    libcurl4-openssl-dev \
    libnlopt-dev \
    libpng-dev \
    git \
    wget \
    libgsl0-dev \
    libgl1-mesa-dev \
    libglu1-mesa-dev \
  && wget https://sourceforge.net/projects/mcmc-jags/files/JAGS/4.x/Source/JAGS-4.2.0.tar.gz -O jags.tar.gz \
  && tar -xf jags.tar.gz \
  && cd JAGS* && ./configure && make && make install \
  && cd / && rm -rf jags.tar.gz JAGS* \
  && R -e 'install.packages(c("devtools")); devtools::install_github("cboettig/pdg_control")'

RUN R -e 'devtools::install(dep=TRUE)'

