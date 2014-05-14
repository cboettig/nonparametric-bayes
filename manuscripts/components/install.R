#!/usr/bin/Rscript

if(!require("nonparametricbayes")){
  install.packages("devtools")
  library("devtools")
  install_github("reshape")
  install_github("cboettig/pdg_control")
  install.packages("..", repos =NULL, dependencies = c("Depends", "Imports", "Suggests"))
}
