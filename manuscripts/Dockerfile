FROM rocker/r-ver:3.2.0
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
#    pandoc \ 
#    pandoc-citeproc \
    libgsl0-dev \
    libgl1-mesa-dev \
    libglu1-mesa-dev
# wget https://sourceforge.net/projects/mcmc-jags/files/JAGS/4.x/Source/JAGS-4.2.0.tar.gz -O jags.tar.gz \
RUN wget https://sourceforge.net/projects/mcmc-jags/files/JAGS/3.x/Source/JAGS-3.4.0.tar.gz/download  -O jags.tar.gz \
  && tar -xf jags.tar.gz \
  && cd JAGS* && ./configure && make && make install \
  && cd / && rm -rf jags.tar.gz JAGS* 
RUN . /etc/environment \
  && R -e 'install.packages(c("devtools")); devtools::install_github("cboettig/pdg_control")' \
  && R -e 'devtools::install_github(c("cboettig/nonparametric-bayes", "cboettig/cboettigR@8bb11be44195e70dfd731cd98f11db1b73e62303"), dep=TRUE)'


## Not needed for code, but only to compile Rmd to pdf:
RUN apt-get update && apt-get -y install \
  texlive-latex-recommended texlive-fonts-recommended 

## debian:jessie version of pandoc-citeproc is too old. 
## also note template is not compatible with newer versions of pandoc. 
RUN apt-get update && apt-get -y install \
  cabal-install libghc-cabal-dev cpphs libghc-pretty-show-dev libexpat1-dev 
RUN cabal update && cabal install hsb2hs
RUN cabal install pandoc-1.13 --force-reinstalls
RUN cabal install pandoc-citeproc-0.5 
RUN ln -s /root/.cabal/bin/pandoc /usr/bin/pandoc
RUN ln -s /root/.cabal/bin/pandoc-citeproc /usr/bin/pandoc-citeproc

COPY . /nonparametric-bayes
WORKDIR /nonparametric-bayes

