
library(knitr)
knit(system.file("examples", "gaussian-process-basics.Rmd", package="nonparametricbayes"))
knit(system.file("examples", "gp_mcmc_example.Rmd", package="nonparametricbayes"))
expect_equal(length(grep("Error", readLines("gp_mcmc_example.md"))), 0)
expect_equal(length(grep("Error", readLines("gaussian-process-basics.md"))), 0)
