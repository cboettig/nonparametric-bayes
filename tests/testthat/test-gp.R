
library(knitr)
knit("gaussian-process-basics.Rmd")
knit("gp_mcmc_example.Rmd")
expect_equals(length(grep("Error", readLines("gp_mcmc_example.md"))), 0)
expect_equals(grep("Error", readLines("gaussian-process-basics.md"))), 0)
