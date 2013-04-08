## @knitr libraries
library(pdgControl)
library(nonparametricbayes)
library(reshape2)
library(ggplot2)
library(data.table)
library(tgp)
library(MCMCpack)
library(plyr)
library(knitcitations)


## @knitr plotting-options
theme_set(theme_bw(base_size=10))
theme_update(panel.background = element_rect(fill = "transparent",colour = NA),
             plot.background = element_rect(fill = "transparent",colour = NA))
cbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")


## @knitr stateeq
f <- RickerAllee
p <- c(2, 10, 5) 
K <- 10
allee <- 5


## @knitr sdp-pars
sigma_g <- 0.05
sigma_m <- 0.0
z_g <- function() rlnorm(1, 0, sigma_g)
z_m <- function() 1+(2*runif(1, 0,  1)-1) * sigma_m
x_grid <- seq(0, 1.5 * K, length=101)
h_grid <- x_grid
profit <- function(x,h) pmin(x, h)
delta <- 0.01
OptTime <- 20  # stationarity with unstable models is tricky thing
reward <- 0
xT <- 0
seed_i <- 1
Xo <- K # observations start from
x0 <- Xo # simulation under policy starts from
Tobs <- 35


## @knitr obs
  obs <- sim_obs(Xo, z_g, f, p, Tobs=Tobs, nz=15, 
                 harvest = sort(rep(seq(0, .5, length=7), 5)), seed = seed_i)


## @knitr mle
alt <- par_est(obs,  init = c(r = p[1], 
                              K = mean(obs$x[obs$x>0]), 
                              s = sigma_g))
est <- par_est_allee(obs, f, p,  
                     init = c(r = p[1] + 1, 
                              K = p[2] + 2, 
                              C = p[3] + 2, 
                              s = sigma_g))


## @knitr gp-priors
s2.p <- c(5,5)  
tau2.p <- c(5,1)
d.p = c(10, 1/0.1, 10, 1/0.1)
nug.p = c(10, 1/0.1, 10, 1/0.1) # gamma mean
s2_prior <- function(x) dinvgamma(x, s2.p[1], s2.p[2])
tau2_prior <- function(x) dinvgamma(x, tau2.p[1], tau2.p[2])
d_prior <- function(x) dgamma(x, d.p[1], scale = d.p[2]) + dgamma(x, d.p[3], scale = d.p[4])
nug_prior <- function(x) dgamma(x, nug.p[1], scale = nug.p[2]) + dgamma(x, nug.p[3], scale = nug.p[4])
beta0_prior <- function(x, tau) dnorm(x, 0, tau)
beta = c(0)
priors <- list(s2 = s2_prior, tau2 = tau2_prior, beta0 = dnorm, nug = nug_prior, d = d_prior, ldetK = function(x) 0)


## @knitr gp
  gp <- bgp(X=obs$x, XX=x_grid, Z=obs$y, verb=0,
          meanfn="constant", bprior="b0", BTE=c(2000,16000,2),
          m0r1=FALSE, corr="exp", trace=TRUE, 
          beta = beta, s2.p = s2.p, d.p = d.p, nug.p = nug.p, tau2.p = tau2.p,
          s2.lam = "fixed", d.lam = "fixed", nug.lam = "fixed", tau2.lam = "fixed")      


## @knitr opt
  OPT <- optimal_policy(gp, f, est$f, alt$f,
                        p, est$p, alt$p,
                        x_grid, h_grid, sigma_g, 
                        sigma_g, sigma_g, # est$sigma_g, alt$sigma_g, but those ests are poor
                        delta, xT, profit, reward, OptTime)


## @knitr sim
dt <- simulate_opt(OPT, f, p, x_grid, h_grid, x0, z_g, profit)


## @knitr gp_plot
tgp_dat <- 
    data.frame(  x = gp$XX[[1]], 
                 y = gp$ZZ.km, 
                 ymin = gp$ZZ.km - 2 * sqrt(gp$ZZ.ks2), 
                 ymax = gp$ZZ.km + 2 * sqrt(gp$ZZ.ks2),
                 ymin2 = gp$ZZ.mean - 2 * sqrt(gp$ZZ.vark), 
                 ymax2 = gp$ZZ.mean + 2 * sqrt(gp$ZZ.vark))
  true <- sapply(x_grid, f, 0, p)
  alt_mean <- sapply(x_grid, alt$f, 0, alt$p)
  est_mean <- sapply(x_grid, est$f, 0, est$p)
  models <- data.frame(x=x_grid, GP=tgp_dat$y, 
                       Parametric=est_mean, 
                       True=true, 
                       Structural=alt_mean)
  models <- melt(models, id="x")
  names(models) <- c("x", "method", "value")

ggplot(tgp_dat) + 
  geom_ribbon(aes(x,y,ymin=ymin,ymax=ymax), fill="gray80") +
  geom_ribbon(aes(x,y,ymin=ymin2,ymax=ymax2), fill="gray60") +
  geom_line(data=models, aes(x, value, col=method), alpha=0.8, lwd=1) + 
  geom_point(data=obs, aes(x,y)) + 
  xlab(expression(X[t])) + ylab(expression(X[t+1])) +
  scale_colour_manual(values=cbPalette) 


## @knitr policies_plot
policies <- 
  melt(data.frame(stock=x_grid, 
                  GP = x_grid[OPT$gp_D], 
                  Parametric = x_grid[OPT$est_D],
                  True = x_grid[OPT$true_D],
                  Structural = x_grid[OPT$alt_D]),
       id="stock")
names(policies) <- c("stock", "method", "value")

ggplot(policies, aes(stock, stock - value, color=method)) +
    geom_line(alpha=0.7, lwd=1) + 
    xlab("stock size") + ylab("escapement")  +
    scale_colour_manual(values=cbPalette)


## @knitr sim_plot
ggplot(dt) + 
    geom_line(aes(time, fishstock, 
                  group=interaction(reps, method), 
                  color=method), alpha=.1) +
    scale_colour_manual(values=cbPalette, 
                        guide = guide_legend(override.aes = list(alpha = 1)))


## @knitr posteriors
hyperparameters <- c("index", "s2", "tau2", "beta0", "nug", "d", "ldetK")
posteriors <- melt(gp$trace$XX[[1]][,hyperparameters], id="index")
  prior_curves <- ddply(posteriors, "variable", function(dd){
  grid <- seq(min(dd$value), max(dd$value), length = 100)
  data.frame(value = grid, density = priors[[dd$variable[1]]](grid))
})

ggplot(posteriors) + 
    geom_histogram(aes(x=value, y=..density..), alpha=0.7) +
    geom_line(data=prior_curves, aes(x=value, y=density), col="red") +
    facet_wrap(~ variable, scale="free")


