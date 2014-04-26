
## ----plot-options, cache=FALSE, message=FALSE, warning=FALSE, include=FALSE----
library(nonparametricbayes) 
opts_chunk$set(external=TRUE)
read_chunk("../examples/BUGS/external-chunks.R")

opts_chunk$set(tidy=FALSE, warning=FALSE, message=FALSE, cache=1, 
               comment=NA, verbose=TRUE, fig.width=4, fig.height=3)
 
# Name the cache path and fig.path based on filename...
opts_chunk$set(fig.path = paste("figure/",
                                gsub(".Rmd", "", knitr:::knit_concord$get('infile')),
                                "-", sep=""),
               cache.path = paste("cache/", 
                                  gsub(".Rmd", "", knitr:::knit_concord$get('infile') ), 
                                "/", sep=""))

toggle = "markup" #"hide"  # toggles `results` chunk option 

library(ggplot2) 
theme_set(theme_bw(base_size=12))
cbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", 
               "#F0E442", "#0072B2", "#D55E00", "#CC79A7")


## ----posterior-mode------------------------------------------------------


## ----stateeq-------------------------------------------------------------


## ----sdp-pars, dependson="stateeq"---------------------------------------


## ----obs, dependson="sdp-pars"-------------------------------------------


## ----mle, dependson="obs"------------------------------------------------


## ----mle-output, dependson="mle", results=toggle-------------------------


## ----gp-priors-----------------------------------------------------------


## ----gp, dependson=c("gp-priors", "obs")---------------------------------


## ----gp_traces_densities, dependson="gp"---------------------------------


## ----gp-output, dependson="gp", results=toggle---------------------------


## ----jags-setup----------------------------------------------------------


## ----common-priors-------------------------------------------------------


## ----allen-model---------------------------------------------------------


## ----allen-priors, dependson="common-priors"-----------------------------


## ----allen_priors_xtable, include=TRUE, echo=FALSE-----------------------
allen_priors_xtable <- data.frame(parameter = c("r0", "K", "theta", "sigma"),
           "lower bound" = c(r0_prior_p[1], K_prior_p[1], theta_prior_p[1], stdQ_prior_p[1]),
           "upper bound" = c(r0_prior_p[2], K_prior_p[2], theta_prior_p[2], stdQ_prior_p[2]))


## ----allen-mcmc, dependson=c("allen-model", "allen-pars", "jags-setup"), results=toggle----


## ----allen-traces, dependson="allen-mcmc"--------------------------------


## ----allen-posteriors, dependson=c("allen-traces", "allen-priors")-------


## ----allen-output, dependson=c("posterior-mode", "allen-traces"), results=toggle----


## ----ricker-model--------------------------------------------------------


## ----ricker-priors, dependson="common-priors"----------------------------


## ----ricker_priors_xtable, include=TRUE, echo=FALSE, cache=FALSE---------
ricker_priors_xtable <- data.frame(
  parameter = c("r0", "K", "sigma"),
  "lower bound" = c(r0_prior_p[1], K_prior_p[1], stdQ_prior_p[1]),
  "upper bound" = c(r0_prior_p[2], K_prior_p[2], stdQ_prior_p[2]))


## ----ricker-mcmc, dependson="ricker-model", results=toggle---------------


## ----ricker-traces, dependson="ricker-mcmc"------------------------------


## ----ricker-posteriors, dependson=c("ricker-traces", "ricker-priors")----


## ----ricker-output, dependson=c("posterior-mode", "ricker-traces"), results=toggle----


## ----myers-model---------------------------------------------------------


## ----myers-priors--------------------------------------------------------


## ----myers_prior_xtable, include=TRUE, echo=FALSE------------------------
myers_priors_xtable <- data.frame(parameter = c("r0", "K", "theta", "sigma"),
           "lower bound" = c(r0_prior_p[1], K_prior_p[1], theta_prior_p[1], stdQ_prior_p[1]),
           "upper bound" = c(r0_prior_p[2], K_prior_p[2], theta_prior_p[2], stdQ_prior_p[2]))


## ----myers-mcmc, dependson="myers-model", results=toggle-----------------


## ----myers-traces, dependson="myers-mcmc"--------------------------------


## ----myers-posteriors, dependson="myers-traces"--------------------------


## ----myers-output, dependson=c("posterior-mode", "myers-traces"), results=toggle----


## ----assemble-models, dependson=c("myers-output", "ricker-output", "allen-output", "gp-output", "mle-output")----


## ----gp-opt, dependson="gp-output"---------------------------------------


## ----mle-opt, dependson="mle-output"-------------------------------------


## ----allen-opt, dependson="allen-output"---------------------------------


## ----ricker-opt, dependson="ricker-output"-------------------------------


## ----myers-opt, dependson="myers-output"---------------------------------


## ----assemble-opt, dependson=c("gp-opt", "mle-opt", "allen-opt", "ricker-opt", "myers-opt")----


## ----sims, dependson="assemble-opt"--------------------------------------


## ----profits, dependson="sims", results=toggle---------------------------


## ----dic_calc, dependson=c("posterior-mode", "myers-output", "ricker-output", "allen-output", "gp-output", "mle-output"), include=FALSE, echo=FALSE----


## ----Figureb_posteriors, dependson=c("assemble-models", "par-fns"), include=TRUE, echo=FALSE, fig.cap="Points show the training data of stock-size over time.  Curves show the posterior step-ahead predictions based on each of the estimated models.", fig.width=8, fig.height=6----


## ----deviances, dependson=c("dic_calc"), include=FALSE, echo=FALSE-------


## ----Table1, dependson="deviances", include=TRUE, results="asis", echo=FALSE, cache=FALSE----
xtable::xtable(rbind(dictable, aictable, bictable), caption="Model choice scores for several common criteria all (wrongly) select the simplest model. As the true (Allen) model is not distinguishable from the simpler (Ricker) model in the region of the observed data, this error cannot be avoided regardless of the model choice criterion. This highlights the danger of model choice when the selected model will be used outside of the observed range of the data.")


## ----statespace_posteriors, dependson=c("assemble-models"), fig.cap="Graph of the inferred Gaussian process compared to the true process and maximum-likelihood estimated process.  Graph shows the expected value for the function $f$ under each model.  Two standard deviations from the estimated Gaussian process covariance with (light grey) and without (darker grey) measurement error are also shown.  The training data is also shown as black points.  The GP is conditioned on (0,0), shown as a pseudo-data point.", include=TRUE, echo=FALSE, fig.width=8, fig.height=6----


## ----out_of_sample_predictions, include=TRUE, echo=FALSE, fig.width=8, fig.height=6, fig.cap="Out of sample predictions of the dynamics under each model.  Points show the stock size simulated by the true model.  Overlay shows the range of states predicted by each model, based on the state observed in the previous time step. The Ricker model always predicts population growth, while the actual population shrinks in each step as the initial condition falls below the Allee threshold of the underlying model (Allen).  Note that the GP is both more pessimistic and more uncertain about the future state than the parameteric models, while the realized state often falls outside of the expected range forecasted by the structurally incorrect Myers and Ricker models."----


## ----Figure2, dependson = c("assemble-opt"), fig.cap="The steady-state optimal policy (infinite boundary) calculated under each model.  Policies are shown in terms of target escapement, $S_t$, as under models such as this a constant escapement policy is expected to be optimal [@Reed1979].", include=TRUE, echo=FALSE----


## ----Figure3, dependson=c("sim"), fig.cap="Gaussian process inference outperforms parametric estimates. Shown are 100 replicate simulations of the stock dynamics (eq 1) under the policies derived from each of the estimated models, as well as the policy based on the exact underlying model.",  include=TRUE, echo=FALSE----


## ----Figure4, dependson=c("profits"), fig.cap="Histograms of the realized net present value of the fishery over a range of simulated data and resulting parameter estimates. For each data set, the three models are estimated as described above. Values plotted are the averages of a given policy over 100 replicate simulations. Details and code provided in the supplement.",  include=TRUE, echo=FALSE, fig.width=6, fig.height=6----

ggplot(actual_over_optimal, aes(value)) + geom_histogram(aes(fill=variable)) + 
  facet_wrap(~variable, ncol=2)  + 
  guides(legend.position = "none") +
  xlab("Total profit by replicate") + 
  scale_fill_manual(values=colorkey) # density plots fail when delta fn

# ggplot(actual_over_optimal, aes(value)) + geom_histogram(aes(fill=variable), binwidth=0.1) + 
#  xlab("Total profit by replicate")+ scale_fill_manual(values=colorkey)
# ggplot(actual_over_optimal, aes(value, fill=variable, color=variable)) + # density plots fail when delta fn
#  stat_density(aes(y=..density..), position="stack", adjust=3, alpha=.9) + 
#  xlab("Total profit by replicate")+ scale_fill_manual(values=colorkey)+ scale_color_manual(values=colorkey)


## ----, echo=FALSE, include=TRUE, fig.height=6, fig.cap="Traces from the MCMC estimate of the Ricker model", fig.width=6----
plot_ricker_traces


## ----, echo=FALSE, include=TRUE, fig.cap="Posteriors from the MCMC estimate of the Ricker model", fig.width=6, fig.height=4----
ggplot(ricker_posteriors, aes(value)) + 
  stat_density(geom="path", position="identity", alpha=0.7) +
  facet_wrap(~ variable, scale="free", ncol=2)


## ----echo=FALSE, results = "asis", include = TRUE, cache=FALSE-----------
xtable::xtable(ricker_priors_xtable,
  caption = "parameterization range for the uniform priors in the Ricker model")


## ----, echo=FALSE, include=TRUE, fig.height=6, fig.cap="Traces from the MCMC estimate of the Myers model", fig.width=6----
plot_myers_traces


## ----, echo=FALSE, include=TRUE, fig.cap="Posterior distributions from the MCMC estimates of the Myers model", fig.width=6, fig.height=6----
ggplot(myers_posteriors, aes(value)) + 
  stat_density(geom="path", position="identity", alpha=0.7) +
  facet_wrap(~ variable, scale="free", ncol=2)


## ----echo=FALSE, results="asis", include = TRUE, cache=FALSE-------------
xtable::xtable(myers_priors_xtable,
           caption = "parameterization range for the uniform priors in the Myers model")


## ----, echo=FALSE, include=TRUE, fig.height=6, fig.cap="Traces from the MCMC estimate of the Allen model", fig.width=6----
plot_allen_traces


## ----, echo=FALSE, include=TRUE, fig.cap="Posteriors from the MCMC estimate of the Allen model", fig.width=6, fig.height=6----
ggplot(allen_posteriors, aes(value)) + 
  stat_density(geom="path", position="identity", alpha=0.7) +
  facet_wrap(~ variable, scale="free", ncol=2)


## ----echo=FALSE, results = "asis", include=TRUE, cache=FALSE-------------
xtable::xtable(allen_priors_xtable,
  caption = "parameterization range for the uniform priors in the Allen model")


## ----echo=FALSE, include=TRUE, fig.cap = "Traces from the MCMC estimates of the GP model"----
gp_assessment_plots$traces_plot


## ----echo=FALSE, include=TRUE, fig.cap="Posterior distributions from the MCMC estimate of the GP model. Prior curves shown in red."----
gp_assessment_plots$posteriors_plot


