
## ----plot-options, message=FALSE, warning=FALSE, include=FALSE, echo=FALSE----
cbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", 
               "#F0E442", "#0072B2", "#D55E00", "#CC79A7")



## ----libraries, cache=FALSE, message=FALSE, warning=FALSE, include=FALSE, echo=FALSE----
library("methods")
library("knitr")

library("nonparametricbayes") 
library("pdgControl")

library("ggplot2") 
library("reshape2")
library("plyr")
library("data.table")

library("R2jags")
library("emdbook") # for as.mcmc.bugs (?)
library("coda")  # for as.mcmc

library("modeest")
library("MASS")
library("pander")

opts_chunk$set(fig.width=4, fig.height=3, echo=FALSE)
theme_set(theme_bw(base_size=12))
toggle = "hide" # results argument


## ----posterior-mode------------------------------------------------------
library(modeest)
posterior.mode <- function(x) {
  mlv(x, method="shorth")$M
}

## ----stateeq-------------------------------------------------------------
f <- RickerAllee
p <- c(2, 8, 2)
K <- 8
allee <- 2 

## ----sdp-pars, dependson="stateeq"---------------------------------------
sigma_g <- 0.05
sigma_m <- 0.0
z_g <- function() rlnorm(1, 0, sigma_g)
z_m <- function() 1
x_grid <- seq(0, 1.5 * K, length=50)
h_grid <- x_grid
profit <- function(x,h) pmin(x, h)
delta <- 0.01
OptTime <- 50  # stationarity with unstable models is tricky thing
reward <- 0
xT <- 0
Xo <-  allee+.5# observations start from
x0 <- 10 # simulation under policy starts from
Tobs <- 40
MaxT = 1000 # timeout for value iteration convergence

## ----obs, dependson="sdp-pars", fig.keep='none'--------------------------
  set.seed(1234)
  #harvest <- sort(rep(seq(0, .5, length=7), 5))
  x <- numeric(Tobs)
  x[1] <- Xo
  nz <- 1
  for(t in 1:(Tobs-1))
    x[t+1] = z_g() * f(x[t], h=0, p=p)
  obs <- data.frame(x = c(rep(0,nz), 
                          pmax(rep(0,Tobs-1), x[1:(Tobs-1)])), 
                    y = c(rep(0,nz), 
                          x[2:Tobs]))
raw_plot <- ggplot(data.frame(time = 1:Tobs, x=x), aes(time,x)) + geom_line()


## ----mle, dependson="obs"------------------------------------------------
set.seed(12345)
estf <- function(p){ 
    mu <- f(obs$x,0,p)
    -sum(dlnorm(obs$y, log(mu), p[4]), log=TRUE)
}
par <- c(p[1]*rlnorm(1,0,.1), 
         p[2]*rlnorm(1,0,.1), 
         p[3]*rlnorm(1,0, .1), 
         sigma_g * rlnorm(1,0,.1))
o <- optim(par, estf, method="L", lower=c(1e-5,1e-5,1e-5,1e-5))
f_alt <- f
p_alt <- c(as.numeric(o$par[1]), as.numeric(o$par[2]), as.numeric(o$par[3]))
sigma_g_alt <- as.numeric(o$par[4])

est <- list(f = f_alt, p = p_alt, sigma_g = sigma_g_alt, mloglik=o$value)


## ----mle-output, dependson="mle", results=toggle-------------------------
true_means <- sapply(x_grid, f, 0, p)
est_means <- sapply(x_grid, est$f, 0, est$p)


## ----gp-priors-----------------------------------------------------------
s2.p <- c(5,5)  
d.p = c(10, 1/0.1)

## ----gp, dependson=c("gp-priors", "obs")---------------------------------
gp <- gp_mcmc(obs$x, y=obs$y, n=1e5, s2.p = s2.p, d.p = d.p)
gp_dat <- gp_predict(gp, x_grid, burnin=1e4, thin=300)

## ----gp_traces_densities, dependson="gp"---------------------------------
gp_assessment_plots <- summary_gp_mcmc(gp, burnin=1e4, thin=300)

## ----gp-output, dependson="gp", results=toggle---------------------------
tgp_dat <- 
    data.frame(  x = x_grid, 
                 y = gp_dat$E_Ef, 
                 ymin = gp_dat$E_Ef - 2 * sqrt(gp_dat$E_Vf), 
                 ymax = gp_dat$E_Ef + 2 * sqrt(gp_dat$E_Vf) )


## ----jags-setup----------------------------------------------------------
y <- x 
N <- length(x);
jags.data <- list("N","y")
n.chains <- 6
n.iter <- 1e6
n.burnin <- floor(10000)
n.thin <- max(1, floor(n.chains * (n.iter - n.burnin)/1000))
n.update <- 10


## ----common-priors-------------------------------------------------------
stdQ_prior_p <- c(1e-6, 100)
stdR_prior_p <- c(1e-6, .1)
stdQ_prior  <- function(x) dunif(x, stdQ_prior_p[1], stdQ_prior_p[2])
stdR_prior  <- function(x) dunif(x, stdR_prior_p[1], stdR_prior_p[2])


## ----allen-model---------------------------------------------------------
K_prior_p <- c(0.01, 20.0)
r0_prior_p <- c(0.01, 6.0)
theta_prior_p <- c(0.01, 20.0)

bugs.model <- 
paste(sprintf(
"model{
  K     ~ dunif(%s, %s)
  r0    ~ dunif(%s, %s)
  theta ~ dunif(%s, %s)
  stdQ ~ dunif(%s, %s)", 
  K_prior_p[1], K_prior_p[2],
  r0_prior_p[1], r0_prior_p[2],
  theta_prior_p[1], theta_prior_p[2],
  stdQ_prior_p[1], stdQ_prior_p[2]),

  "
  iQ <- 1 / (stdQ * stdQ);
  y[1] ~ dunif(0, 10)
  for(t in 1:(N-1)){
    mu[t] <- log(y[t]) + r0 * (1 - y[t]/K)* (y[t] - theta) / K 
    y[t+1] ~ dlnorm(mu[t], iQ) 
  }
}")
writeLines(bugs.model, "allen_process.bugs")


## ----allen-priors, dependson="common-priors"-----------------------------
K_prior     <- function(x) dunif(x, K_prior_p[1], K_prior_p[2])
r0_prior <- function(x) dunif(x, r0_prior_p[1], r0_prior_p[2])
theta_prior <- function(x) dunif(x, theta_prior_p[1], theta_prior_p[2])
par_priors  <- list(K = K_prior, deviance = function(x) 0 * x, 
                    r0 = r0_prior, theta = theta_prior,
                    stdQ = stdQ_prior)



## ----allen_priors_xtable, include=TRUE, echo=FALSE-----------------------
allen_priors_xtable <- data.frame(parameter = c("$r$", "$K$", "$X_C$", "$\\sigma$"),
           "lower bound" = c(r0_prior_p[1], K_prior_p[1], theta_prior_p[1], stdQ_prior_p[1]),
           "upper bound" = c(r0_prior_p[2], K_prior_p[2], theta_prior_p[2], stdQ_prior_p[2]))


## ----allen-mcmc, dependson=c("allen-model", "allen-pars", "jags-setup"), results=toggle----
jags.params=c("K","r0","theta","stdQ") # be sensible about the order here
jags.inits <- function(){
  list("K"= 10 * rlnorm(1,0, 0.1),
       "r0"= 1 * rlnorm(1,0, 0.1) ,
       "theta"=   5 * rlnorm(1,0, 0.1) , 
       "stdQ"= abs( 0.1 * rlnorm(1,0, 0.1)),
       .RNG.name="base::Wichmann-Hill", .RNG.seed=123)
}

set.seed(1234)
# parallel refuses to take variables as arguments (e.g. n.iter = 1e5 works, but n.iter = n doesn't)
allen_jags <- do.call(jags, list(data=jags.data, inits=jags.inits, 
                                      jags.params, n.chains=n.chains, 
                                      n.iter=n.iter, n.thin=n.thin, 
                                      n.burnin=n.burnin, 
                                      model.file="allen_process.bugs"))

# Run again iteratively if we haven't met the Gelman-Rubin convergence criterion
recompile(allen_jags) # required for parallel
allen_jags <- do.call(autojags, 
											list(object=allen_jags, n.update=n.update, 
                           n.iter=n.iter, n.thin = n.thin))



## ----allen-traces, dependson="allen-mcmc"--------------------------------
tmp <- lapply(as.mcmc(allen_jags), as.matrix) # strip classes the hard way...
allen_posteriors <- melt(tmp, id = colnames(tmp[[1]])) 
names(allen_posteriors) = c("index", "variable", "value", "chain")
plot_allen_traces <- ggplot(allen_posteriors) + geom_line(aes(index, value)) + 
  facet_wrap(~ variable, scale="free", ncol=1)


## ----allen-posteriors, dependson=c("allen-traces", "allen-priors")-------
allen_priors <- ddply(allen_posteriors, "variable", function(dd){
    grid <- seq(min(dd$value), max(dd$value), length = 100) 
    data.frame(value = grid, density = par_priors[[dd$variable[1]]](grid))
})
plot_allen_posteriors <- ggplot(allen_posteriors, aes(value)) + 
  stat_density(geom="path", position="identity", alpha=0.7) +
#  geom_line(data=allen_priors, aes(x=value, y=density), col="red") +  
  facet_wrap(~ variable, scale="free", ncol=3)


## ----allen-output, dependson=c("posterior-mode", "allen-traces"), results=toggle----
A <- allen_posteriors
A$index <- A$index + A$chain * max(A$index) # Combine samples across chains by renumbering index 
pardist <- acast(A, index ~ variable)
bayes_coef <- apply(pardist,2, posterior.mode) 
bayes_pars <- unname(c(bayes_coef["r0"], bayes_coef["K"], bayes_coef["theta"])) # parameters formatted for f
allen_f <- function(x,h,p) unname(RickerAllee(x,h, unname(p[c("r0", "K", "theta")])))
allen_means <- sapply(x_grid, f, 0, bayes_pars)
bayes_pars
head(pardist)



## ----ricker-model--------------------------------------------------------
K_prior_p <- c(0.01, 40.0)
r0_prior_p <- c(0.01, 20.0)
bugs.model <- 
paste(sprintf(
"model{
  K    ~ dunif(%s, %s)
  r0    ~ dunif(%s, %s)
  stdQ ~ dunif(%s, %s)", 
  K_prior_p[1], K_prior_p[2],
  r0_prior_p[1], r0_prior_p[2],
  stdQ_prior_p[1], stdQ_prior_p[2]),
  "
  iQ <- 1 / (stdQ * stdQ);
  y[1] ~ dunif(0, 10)
  for(t in 1:(N-1)){
    mu[t] <- log(y[t]) + r0 * (1 - y[t]/K) 
    y[t+1] ~ dlnorm(mu[t], iQ) 
  }
}")
writeLines(bugs.model, "ricker_process.bugs")


## ----ricker-priors, dependson="common-priors"----------------------------
K_prior     <- function(x) dunif(x, K_prior_p[1], K_prior_p[2])
r0_prior <- function(x) dunif(x, r0_prior_p[1], r0_prior_p[2])
par_priors <- list(K = K_prior, deviance = function(x) 0 * x, 
                   r0 = r0_prior, stdQ = stdQ_prior)



## ----ricker_priors_xtable, include=TRUE, echo=FALSE, cache=FALSE---------
ricker_priors_xtable <- data.frame(
  parameter = c("$r$", "$K$", "$\\sigma$"),
  "lower bound" = c(r0_prior_p[1], K_prior_p[1], stdQ_prior_p[1]),
  "upper bound" = c(r0_prior_p[2], K_prior_p[2], stdQ_prior_p[2]))


## ----ricker-mcmc, dependson="ricker-model", results=toggle---------------
jags.params=c("K","r0", "stdQ")
jags.inits <- function(){
  list("K"= 10 * rlnorm(1,0,.5),
       "r0"= rlnorm(1,0,.5),
       "stdQ"=sqrt(0.05) * rlnorm(1,0,.5),
       .RNG.name="base::Wichmann-Hill", .RNG.seed=123)
}
set.seed(12345) 
ricker_jags <- do.call(jags, 
                       list(data=jags.data, inits=jags.inits, 
                            jags.params, n.chains=n.chains, 
                            n.iter=n.iter, n.thin=n.thin, n.burnin=n.burnin,
                            model.file="ricker_process.bugs"))
recompile(ricker_jags)
ricker_jags <- do.call(autojags, 
                       list(object=ricker_jags, n.update=n.update, 
														n.iter=n.iter, n.thin = n.thin, 
														progress.bar="none"))


## ----ricker-traces, dependson="ricker-mcmc"------------------------------
tmp <- lapply(as.mcmc(ricker_jags), as.matrix) # strip classes the hard way...
ricker_posteriors <- melt(tmp, id = colnames(tmp[[1]])) 
names(ricker_posteriors) = c("index", "variable", "value", "chain")
plot_ricker_traces <- ggplot(ricker_posteriors) + geom_line(aes(index, value)) + 
  facet_wrap(~ variable, scale="free", ncol=1)

## ----ricker-posteriors, dependson=c("ricker-traces", "ricker-priors")----
ricker_priors <- ddply(ricker_posteriors, "variable", function(dd){
    grid <- seq(min(dd$value), max(dd$value), length = 100) 
    data.frame(value = grid, density = par_priors[[dd$variable[1]]](grid))
})
# plot posterior distributions
plot_ricker_posteriors <- ggplot(ricker_posteriors, aes(value)) + 
  stat_density(geom="path", position="identity", alpha=0.7) +
#  geom_line(data=ricker_priors, aes(x=value, y=density), col="red") +  # don't plot priors 
  facet_wrap(~ variable, scale="free", ncol=2)

## ----ricker-output, dependson=c("posterior-mode", "ricker-traces"), results=toggle----
A <- ricker_posteriors
A$index <- A$index + A$chain * max(A$index) # Combine samples across chains by renumbering index 
ricker_pardist <- acast(A, index ~ variable)
bayes_coef <- apply(ricker_pardist,2, posterior.mode) 
ricker_bayes_pars <- unname(c(bayes_coef["r0"], bayes_coef["K"]))
ricker_f <- function(x,h,p){
  sapply(x, function(x){ 
    x <- pmax(0, x-h) 
    pmax(0, x * exp(p["r0"] * (1 - x / p["K"] )) )
  })
}
ricker_means <- sapply(x_grid, Ricker, 0, ricker_bayes_pars[c(1,2)])
head(ricker_pardist)
ricker_bayes_pars



## ----myers-model---------------------------------------------------------
r0_prior_p <- c(.0001, 10.0)
theta_prior_p <- c(.0001, 10.0)
K_prior_p <- c(.0001, 40.0)
bugs.model <- 
paste(sprintf(
"model{
  r0    ~ dunif(%s, %s)
  theta    ~ dunif(%s, %s)
  K    ~ dunif(%s, %s)
  stdQ ~ dunif(%s, %s)", 
  r0_prior_p[1], r0_prior_p[2],
  theta_prior_p[1], theta_prior_p[2],
  K_prior_p[1], K_prior_p[2],
  stdQ_prior_p[1], stdQ_prior_p[2]),

  "
  iQ <- 1 / (stdQ * stdQ);

  y[1] ~ dunif(0, 10)
  for(t in 1:(N-1)){
    mu[t] <- log(r0)  + theta * log(y[t]) - log(1 + pow(abs(y[t]), theta) / K)
    y[t+1] ~ dlnorm(mu[t], iQ) 
  }
}")
writeLines(bugs.model, "myers_process.bugs")


## ----myers-priors--------------------------------------------------------
K_prior     <- function(x) dunif(x, K_prior_p[1], K_prior_p[2])
r_prior     <- function(x) dunif(x, r0_prior_p[1], r0_prior_p[2])
theta_prior <- function(x) dunif(x, theta_prior_p[1], theta_prior_p[2])
par_priors <- list( deviance = function(x) 0 * x, K = K_prior,
                    r0 = r_prior, theta = theta_prior, 
                    stdQ = stdQ_prior)



## ----myers_prior_xtable, include=TRUE, echo=FALSE------------------------
myers_priors_xtable <- data.frame(parameter = c("$r$", "$K$", "$\\theta$", "$\\sigma$"),
           "lower bound" = c(r0_prior_p[1], K_prior_p[1], theta_prior_p[1], stdQ_prior_p[1]),
           "upper bound" = c(r0_prior_p[2], K_prior_p[2], theta_prior_p[2], stdQ_prior_p[2]))


## ----myers-mcmc, dependson="myers-model", results=toggle-----------------
jags.params=c("r0", "theta", "K", "stdQ")
jags.inits <- function(){
  list("r0"= 1 * rlnorm(1,0,.1), 
       "K"=    10 * rlnorm(1,0,.1),
       "theta" = 1 * rlnorm(1,0,.1),  
       "stdQ"= sqrt(0.2) * rlnorm(1,0,.1),
       .RNG.name="base::Wichmann-Hill", .RNG.seed=123)
}
set.seed(12345)
myers_jags <- do.call(jags, 
                      list(data=jags.data, inits=jags.inits, 
													 jags.params, n.chains=n.chains, 
													 n.iter=n.iter, n.thin=n.thin,
                           n.burnin=n.burnin, 
                           model.file="myers_process.bugs"))
recompile(myers_jags)
myers_jags <- do.call(autojags, 
                      list(myers_jags, n.update=n.update, 
                           n.iter=n.iter, n.thin = n.thin, 
                           progress.bar="none"))


## ----myers-traces, dependson="myers-mcmc"--------------------------------
tmp <- lapply(as.mcmc(myers_jags), as.matrix) # strip classes
myers_posteriors <- melt(tmp, id = colnames(tmp[[1]])) 
names(myers_posteriors) = c("index", "variable", "value", "chain")
plot_myers_traces <- ggplot(myers_posteriors) + 
  geom_line(aes(index, value)) + # priors, need to fix order though
  facet_wrap(~ variable, scale="free", ncol=1)


## ----myers-posteriors, dependson="myers-traces"--------------------------
par_prior_curves <- ddply(myers_posteriors, "variable", function(dd){
    grid <- seq(min(dd$value), max(dd$value), length = 100) 
    data.frame(value = grid, density = par_priors[[dd$variable[1]]](grid))
})
plot_myers_posteriors <- ggplot(myers_posteriors, aes(value)) + 
  stat_density(geom="path", position="identity", alpha=0.7) +
#  geom_line(data=par_prior_curves, aes(x=value, y=density), col="red") +  # Whoops, these are misaligned. see table instead 
  facet_wrap(~ variable, scale="free", ncol=3)


## ----myers-output, dependson=c("posterior-mode", "myers-traces"), results=toggle----
A <- myers_posteriors
A$index <- A$index + A$chain * max(A$index) # Combine samples across chains by renumbering index 
myers_pardist <- acast(A, index ~ variable)
bayes_coef <- apply(myers_pardist,2, posterior.mode) # much better estimates
myers_bayes_pars <- unname(c(bayes_coef["r0"], bayes_coef["theta"], bayes_coef["K"]))
myers_means <- sapply(x_grid, Myer_harvest, 0, myers_bayes_pars)
myers_f <- function(x,h,p) Myer_harvest(x, h, p[c("r0", "theta", "K")])
head(myers_pardist)
myers_bayes_pars



## ----assemble-models, dependson=c("myers-output", "ricker-output", "allen-output", "gp-output", "mle-output")----
models <- data.frame(x=x_grid, 
										 GP=tgp_dat$y, 
										 True=true_means, 
                     MLE=est_means, 
										 Ricker=ricker_means, 
                     Allen = allen_means,
                     Myers = myers_means)
models <- melt(models, id="x")
# some labels
names(models) <- c("x", "method", "value")
# labels for the colorkey too
model_names = c("GP", "True", "MLE", "Ricker", "Allen", "Myers")
colorkey=cbPalette
names(colorkey) = model_names 
step_ahead <- function(x, f, p){
  h = 0
  x_predict <- sapply(x, f, h, p)
  n <- length(x_predict) - 1
  y <- c(x[1], x_predict[1:n])
  y
}

## ----gp-opt, dependson="gp-output"---------------------------------------
matrices_gp <- gp_transition_matrix(gp_dat$Ef_posterior, gp_dat$Vf_posterior, x_grid, h_grid) 
opt_gp <- value_iteration(matrices_gp, x_grid, h_grid, MaxT, xT, profit, delta, reward)

## ----mle-opt, dependson="mle-output"-------------------------------------
matrices_true <- f_transition_matrix(f, p, x_grid, h_grid, sigma_g)
opt_true <- value_iteration(matrices_true, x_grid, h_grid, OptTime=MaxT, xT, profit, delta=delta)
matrices_estimated <- f_transition_matrix(est$f, est$p, x_grid, h_grid, est$sigma_g)
opt_estimated <- value_iteration(matrices_estimated, x_grid, h_grid, OptTime=MaxT, xT, profit, delta=delta)


## ----allen-opt, dependson="allen-output"---------------------------------
matrices_allen <- parameter_uncertainty_SDP(allen_f, x_grid, h_grid, pardist, 4)
opt_allen <- value_iteration(matrices_allen, x_grid, h_grid, OptTime=MaxT, xT, profit, delta=delta)

## ----ricker-opt, dependson="ricker-output"-------------------------------
matrices_ricker <- parameter_uncertainty_SDP(ricker_f, x_grid, h_grid, as.matrix(ricker_pardist), 3)
opt_ricker <- value_iteration(matrices_ricker, x_grid, h_grid, OptTime=MaxT, xT, profit, delta=delta)

## ----myers-opt, dependson="myers-output"---------------------------------
matrices_myers <- parameter_uncertainty_SDP(myers_f, x_grid, h_grid, as.matrix(myers_pardist), 4)
myers_alt <- value_iteration(matrices_myers, x_grid, h_grid, OptTime=MaxT, xT, profit, delta=delta)


## ----assemble-opt, dependson=c("gp-opt", "mle-opt", "allen-opt", "ricker-opt", "myers-opt")----
OPT = data.frame(GP = opt_gp$D, True = opt_true$D, MLE = opt_estimated$D, Ricker = opt_ricker$D, Allen = opt_allen$D, Myers = myers_alt$D)
colorkey=cbPalette
names(colorkey) = names(OPT) 

## ----sims, dependson="assemble-opt"--------------------------------------
sims <- lapply(OPT, function(D){
  set.seed(1)
  lapply(1:100, function(i) 
    ForwardSimulate(f, p, x_grid, h_grid, x0, D, z_g, profit=profit, OptTime=OptTime)
  )
})
# turn the list into a data.frame
keep.as.columns <- names(sims[[1]][[1]])
dat <- melt(sims, id=keep.as.columns)
sims_data <- data.table(dat)
setnames(sims_data, c("L1", "L2"), c("method", "reps")) 
# Legend in original ordering please, not alphabetical: 
sims_data$method = factor(sims_data$method, ordered=TRUE, levels=names(OPT))


## ----profits, dependson="sims", results=toggle---------------------------
Profit <- sims_data[, sum(profit), by=c("reps", "method")]
tmp <- dcast(Profit, reps ~ method)
#tmp$Allen <- tmp[,"Allen"] + rnorm(dim(tmp)[1], 0, 1) # jitter for plotting
tmp <- tmp / tmp[,"True"]
tmp <- melt(tmp[2:dim(tmp)[2]])
actual_over_optimal <-subset(tmp, variable != "True")


## ----dic_calc, dependson=c("posterior-mode", "myers-output", "ricker-output", "allen-output", "gp-output", "mle-output"), include=FALSE, echo=FALSE----
dic.dic <- function (x) sum(x$deviance) +  sum(x[[2]])
recompile(allen_jags)
allen_dic <- dic.dic(dic.samples(allen_jags$model, n.iter=1000, type="popt"))
recompile(ricker_jags)
ricker_dic <- dic.dic(dic.samples(ricker_jags$model, n.iter=1000, type="popt"))
recompile(myers_jags)
myers_dic <- dic.dic(dic.samples(myers_jags$model, n.iter=1000, type="popt"))
dictable <- data.frame(Allen = allen_dic + 2*length(bayes_pars), 
                       Ricker = ricker_dic + 2*length(ricker_bayes_pars),
                       Myers = myers_dic + 2*length(myers_bayes_pars), 
                       row.names = c("DIC"))




## ----figure_1, dependson=c("assemble-models", "par-fns", "plot-options"), fig.cap="Points show the training data of stock-size over time.  Curves show the posterior step-ahead predictions based on each of the estimated models. Observe that all models are fitting the data reasonably well.", fig.width=8, fig.height=6----
step_ahead_posteriors <- function(x){
  gp_f_at_obs <- gp_predict(gp, x, burnin=1e4, thin=300)
  df_post <- melt(lapply(sample(100, 30),  
  function(i){
    data.frame(time = 1:length(x), stock = x, 
                GP = mvrnorm(1, gp_f_at_obs$Ef_posterior[,i], gp_f_at_obs$Cf_posterior[[i]]),
                True = step_ahead(x,f,p),  
                MLE = step_ahead(x,f,est$p), 
                Allen = step_ahead(x, allen_f, pardist[i,]), 
                Ricker = step_ahead(x, ricker_f, ricker_pardist[i,]), 
                Myers = step_ahead(x, myers_f, myers_pardist[i,]))
  }), id=c("time", "stock"))
}
df_post <- step_ahead_posteriors(x)
figure1b_posteriors <- ggplot(df_post) + geom_point(aes(time, stock)) + 
  geom_line(aes(time, value, col=variable, group=interaction(L1,variable)), alpha=.1) + 
  facet_wrap(~variable) + 
  scale_colour_manual(values=colorkey, guide = guide_legend(override.aes = list(alpha = 1))) +  
  theme(legend.position="none")
figure1b_posteriors



## ----deviances, dependson=c("dic_calc"), include=FALSE, echo=FALSE-------
allen_deviance  <- - posterior.mode(pardist[,'deviance'])
ricker_deviance <- - posterior.mode(ricker_pardist[,'deviance'])
myers_deviance  <- - posterior.mode(myers_pardist[,'deviance'])
true_deviance   <- 2*estf(c(p, sigma_g))
mle_deviance    <- 2*estf(c(est$p, est$sigma_g))
aictable <- data.frame(Allen = allen_deviance + 2*(1+length(bayes_pars)),  # +1 for noise parameter
                       Ricker = ricker_deviance + 2*(1+length(ricker_bayes_pars)),
                       Myers = myers_deviance + 2*(1+length(myers_bayes_pars)), 
                       row.names = c("AIC"))
bictable <- data.frame(Allen = allen_deviance + log(length(x))*(1+length(bayes_pars)), 
                       Ricker = ricker_deviance + log(length(x))*(1+length(ricker_bayes_pars)),
                       Myers = myers_deviance + log(length(x))*(1+length(myers_bayes_pars)), 
                       row.names = c("BIC"))
xtable::xtable(rbind(dictable, aictable, bictable))


## ----Table1, dependson="deviances", include=TRUE, results="asis", echo=FALSE, cache=FALSE----
xtable::xtable(rbind(dictable, aictable, bictable), caption="Model selection scores for several common criteria all (wrongly) select the simplest model. As the true (Allen) model is not distinguishable from the simpler (Ricker) model in the region of the observed data, this error cannot be avoided regardless of the model choice criterion. This highlights the danger of model choice when the selected model will be used outside of the observed range of the data.")


## ----figure_2, dependson=c("assemble-models", "plot-options"), fig.cap="Graph of the inferred Gaussian process compared to the true process and maximum-likelihood estimated process.  Graph shows the expected value for the function $f$ under each model.  Two standard deviations from the estimated Gaussian process covariance with (light grey) and without (darker grey) measurement error are also shown.  The training data is also shown as black points.  The GP is conditioned on (0,0), shown as a pseudo-data point.", fig.width=8, fig.height=6----
x_grid_short <- x_grid[1:40]
gp_short <- gp_predict(gp, x_grid_short, burnin=1e4, thin=300)
models_posteriors <- 
  melt(lapply(sample(100, 50), 
              function(i){
    sample_gp <- mvrnorm(1, 
                            gp_short$Ef_posterior[,i],         
                            gp_short$Cf_posterior[[i]])
    data.frame(stock = x_grid_short, 
               GP = sample_gp,
               y = sample_gp,
               ymin = sample_gp - 2 * sqrt(gp_short$E_Vf), 
               ymax = sample_gp + 2 * sqrt(gp_short$E_Vf), 
               True = sapply(x_grid_short,f,0, p),  
               MLE = sapply(x_grid_short,f,0, est$p), 
               Allen = sapply(x_grid_short, allen_f, 0, pardist[i,]), 
               Ricker = sapply(x_grid_short, ricker_f, 0, ricker_pardist[i,]), 
               Myers = sapply(x_grid_short, myers_f, 0, myers_pardist[i,]))
             }), 
       id=c("stock", "y", "ymin", "ymax"))
ggplot(models_posteriors) + 
    geom_ribbon(aes(x=stock, y=y, ymin=ymin, ymax=ymax, group=L1), 
                  fill = "gray80", 
                  data=subset(models_posteriors, variable == "GP")) + 
    geom_line(aes(stock, value, col = variable, 
                  group=interaction(L1,variable)), 
              alpha=.2) + 
    geom_point(data = obs, aes(x,y), alpha = 0.8) + 
    xlab(expression(X[t])) + ylab(expression(X[t+1])) +
    facet_wrap(~variable) + 
    scale_colour_manual(values=colorkey) +  
    theme(legend.position="none")


## ----figure_3, fig.width=8, fig.height=6, fig.cap="Beyond the tipping point (outside the range of the training data), the true dynamics fall outside the uncertainty bands of the wrong parametric models, but inside the uncertainty predicted by the GP. Points show the stock size simulated by the true model.  Overlay shows the range of states predicted by each model, based on the state observed in the previous time step. The Ricker model always (wrongly) predicts positive population growth, while the actual population shrinks in each step as the initial condition falls below the Allee threshold of the underlying model (Allen).  Note that the GP is both more pessimistic and more uncertain about the future state than the parametric models.", dependson="plot-options"----
y <- numeric(8)
y[1] <- 4.5
for(t in 1:(length(y)-1))
      y[t+1] = z_g() * f(y[t], h=0, p=p)
# predicts means, does not reflect uncertainty estimate!
crash_data <- step_ahead_posteriors(y)
crash_data <- subset(crash_data, variable %in% c("GP", "Allen", "Ricker", "Myers"))
ggplot(crash_data) + 
  geom_boxplot(aes(as.factor(as.integer(time)), value, 
                   fill = variable, col=variable), 
               alpha=.7, outlier.size=1, position="identity") + 
#  geom_line(aes(time, value, col = variable, 
#            group=interaction(L1,variable)), alpha=.1) + 
  geom_point(aes(time, stock), size = 3) + 
  scale_fill_manual(values=colorkey[c("GP", "Allen", "Ricker", "Myers")], 
                      guide = guide_legend(override.aes = list(alpha = 1))) +  
  scale_colour_manual(values=colorkey[c("GP", "Allen", "Ricker", "Myers")], 
                      guide = guide_legend(override.aes = list(alpha = 1))) +  
  facet_wrap(~variable) + 
  theme(legend.position="none") + xlab("time") + ylab("stock size") 


## ----figure_4, dependson = c("assemble-opt"), fig.cap="The steady-state optimal policy (infinite boundary) calculated under each model.  Policies are shown in terms of target escapement, $S_t$, as under models such as this a constant escapement policy is expected to be optimal [@Reed1979]."----
policies <- melt(data.frame(stock=x_grid, sapply(OPT, function(x) x_grid[x])), id="stock")
names(policies) <- c("stock", "method", "value")

ggplot(policies, aes(stock, stock - value, color=method)) +
  geom_line(lwd=1.2, alpha=0.8) + xlab("stock size") + ylab("escapement")  +
  scale_colour_manual(values=colorkey)


## ----figure_5, dependson=c("sim"), fig.cap="Gaussian process inference outperforms parametric estimates. Shown are 100 replicate simulations of the stock dynamics (eq 1) under the policies derived from each of the estimated models, as well as the policy based on the exact underlying model.", fig.width=6, fig.height=6----
ggplot(sims_data) + 
  geom_line(aes(time, fishstock, group=interaction(reps,method), color=method), alpha=.1) +
  scale_colour_manual(values=colorkey, guide = guide_legend(override.aes = list(alpha = 1))) + 
  facet_wrap(~method) + guides(legend.position="none")



## ----figure_6, dependson=c("profits"), fig.cap="Histograms of the realized net present value of the fishery over a range of simulated data and resulting parameter estimates. For each data set, the three models are estimated as described above. Values plotted are the averages of a given policy over 100 replicate simulations. Details and code provided in the supplement.", fig.width=6, fig.height=6----

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




## ----figure_S1, fig.cap = "Traces from the MCMC estimates of the GP model show reasonable mixing (no trend) and sampling rejection rate (no piecewise jumps)", dependson="plot-options"----
gp_assessment_plots$traces_plot


## ----figure_S2, fig.cap="Posterior distributions from the MCMC estimate of the GP model. Prior curves shown in red."----
gp_assessment_plots$posteriors_plot


## ----figure_S3, fig.height=6, fig.cap="Traces from the MCMC estimates of the Ricker model show reasonable mixing (no trend) and sampling rejection rate (no piecewise jumps)", fig.width=6----
plot_ricker_traces


## ----figure_S4, fig.cap="Posteriors from the MCMC estimate of the Ricker model", fig.width=6, fig.height=4----
ggplot(ricker_posteriors, aes(value)) + 
  stat_density(geom="path", position="identity", alpha=0.7) +
  facet_wrap(~ variable, scale="free", ncol=2)


## ----Table S1, results = "asis"------------------------------------------
pander::pandoc.table(ricker_priors_xtable,
  caption = "Parameterization range for the uniform priors in the Ricker model")


## ----figure_S5, fig.height=6, fig.cap="Traces from the MCMC estimates of the Myers model show reasonable mixing (no trend) and sampling rejection rate (no piecewise jumps)", fig.width=6----
plot_myers_traces


## ----figure_S6, fig.cap="Posterior distributions from the MCMC estimates of the Myers model", fig.width=6, fig.height=6----
ggplot(myers_posteriors, aes(value)) + 
  stat_density(geom="path", position="identity", alpha=0.7) +
  facet_wrap(~ variable, scale="free", ncol=2)


## ----TableS2, results="asis"---------------------------------------------
pander::pandoc.table(myers_priors_xtable,
           caption = "Parameterization range for the uniform priors in the Myers model")


## ----figure_S7, fig.height=6, fig.cap="Traces from the MCMC estimates of the Allen model show reasonable mixing (no trend) and sampling rejection rate (no piecewise jumps)", fig.width=6----
plot_allen_traces


## ----figure_S8, fig.cap="Posteriors from the MCMC estimate of the Allen model", fig.width=6, fig.height=6----
ggplot(allen_posteriors, aes(value)) + 
  stat_density(geom="path", position="identity", alpha=0.7) +
  facet_wrap(~ variable, scale="free", ncol=2)


## ----TableS3, results = "asis"-------------------------------------------
pander::pandoc.table(allen_priors_xtable,
  caption = "Parameterization range for the uniform priors in the Allen model")


## ----include=FALSE-------------------------------------------------------
unlink("ricker_process.bugs")
unlink("allen_process.bugs")
unlink("myers_process.bugs")


save(list=ls(), file="small-allee.rda")
