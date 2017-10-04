FROM rocker/r-ver:3.2.0
MAINTAINER Carl Boettiger cboettig@ropensci.org

## Install software dependencies
RUN apt-get update \
  && apt-get install -y \
    libssl-dev libxml2-dev libcurl4-openssl-dev git \
  && R -e 'install.packages(c("devtools")); devtools::install_github("cboettig/pdg_control")' \
  && R -e 'devtools::install(dep=TRUE)'

