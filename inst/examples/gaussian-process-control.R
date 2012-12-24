
## @knitr libraries
require(pdgControl)
require(nonparametricbayes)
require(reshape2)
require(ggplot2)
require(data.table)
require(tgp)
library(kernlab)
require(MCMCpack)
require(plyr)



## @knitr graphing-options
opts_knit$set(upload.fun = function(file) socialR::notebook.url(file, cp=TRUE, sync=FALSE))
#opts_chunk$set(dev = 'Cairo_pdf', dev.args=list(""))
opts_chunk$set(dev="png", dev.args=list(bg="transparent"))
opts_chunk$set(comment=NA, tidy=FALSE)
theme_set(theme_bw(base_size=16))
theme_update(panel.background = element_rect(fill = "transparent",colour = NA),
             plot.background = element_rect(fill = "transparent",colour = NA))
cbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
#sha <- gsub("^commit ", "", system("git log -n 1", intern=TRUE)[1])
#short_sha <- gsub("(^.{10}).*", "\\1", sha)
#date <- format(Sys.time(), "%Y-%m-%d-%H-%M-%S")
#opts_chunk$set(fig.path = paste("figure/", date, "-", short_sha, "-", sep=""))






## @knitr BevHolt
f <- BevHolt
p <- c(1.5,.05)
K <- (p[1]-1)/p[2] 



## @knitr Myer
f <- Myer_harvest
pars <- c(1, 2, 4.5) 
p <- pars # shorthand 
K <- p[1] * p[3] / 2 + sqrt( (p[1] * p[3]) ^ 2 - 4 * p[3] ) / 2
allee <- p[1] * p[3] / 2 - sqrt( (p[1] * p[3]) ^ 2 - 4 * p[3] ) / 2 # allee threshold
e_star <- (p[1] * sqrt(p[3]) - 2) / 2 ## Bifurcation point 


## @knitr Ricker
f <- Ricker
p <- c(2.5, 10) 
K <- 12.5


## @knitr May
f <- May
p <- c(r = .75, k = 10, a=1.3, H=1, Q = 3)
K <- 8 

## @knitr plot-may
birth <- function(x) p["r"] * (1-  x / p["k"])
death <- function(x) p["a"] * x ^ (p["Q"] - 1) / (x ^ p["Q"] + p["H"])
df <- data.frame(x = x_grid, b = sapply(x_grid, birth), d = sapply(x_grid, death))
ggplot(df) + geom_line(aes(x, b), col = "blue") + geom_line(aes(x,d), col = "red")


## @knitr RickerAllee
f <- RickerAllee
p <- c(r = 3, K = 10, C=4)
K <- p[2]
allee <- p[3]



## @knitr sdp-pars
sigma_g <- 0.05
sigma_m <- 0.2
z_g <- function() rlnorm(1, 0, sigma_g) #1+(2*runif(1, 0,  1)-1)*sigma_g #
x_grid <- seq(0, 1.5 * K, length=101)
h_grid <- x_grid
profit = function(x,h) pmin(x, h)
delta <- 0.01
OptTime = 50 # should be long enough that profit is greater sustainably harvesting
reward = profit(x_grid[length(x_grid)], x_grid[length(x_grid)]) + 1 / (1 - delta) ^ OptTime 
xT <- 0


## @knitr inits
# x_0_observed is starting condition for simulation of the observed data.  
# It should be in preferred state for bistable model, 
# above Allee threshold for Allee model, 
# and near zero for BH or Ricker models
xT <- allee
x_0_observed <- K



## @knitr sim-obs
Tobs <- 30
x <- numeric(Tobs)
x[1] <- x_0_observed
for(t in 1:(Tobs-1))
  x[t+1] = z_g() * f(x[t], h=0, p=p)
plot(x)

## @knitr lag-data
obs <- data.frame(x=c(0,x[1:(Tobs-1)]),y=c(0,x[2:Tobs]))



## @knitr sim-with-harvest
nz <- 10 # weight on 0,0 pt
Tobs <- 40
harvest <- sort(rep(seq(1, 2, length=5), 8))
x <- numeric(Tobs)
x[1] <- x_0_observed
for(t in 1:(Tobs-1))
  x[t+1] = z_g() * f(x[t], h=2, p=p)
plot(x)


obs <- data.frame(x=c(rep(0,nz), pmax(rep(0,Tobs-1),x[1:(Tobs-1)]-harvest[1:Tobs-1])),y=c(rep(0,nz),x[2:Tobs]))
plot(obs$x, obs$y)


## @knitr MLE-RickerAllee
estf <- function(p){ 
  mu <- f(obs$x,0,p)
  -sum(dlnorm(obs$y, log(mu), p["s"]), log=TRUE)
}
par = c(r = p[1] - 1, 
        C = p[2] - 1, 
        K = p[3] + 2, 
        s = sigma_g + abs(rnorm(1,0,.1)))
o <- optim(par, estf, method="L", lower=c(1e-3,1e-3,1e-3, 1e-3))
f_alt <- f
p_alt <- c(as.numeric(o$par[1]), as.numeric(o$par[2]), as.numeric(o$par[3]))
sigma_g_alt <- as.numeric(o$par[4])
p_alt
sigma_g_alt





## @knitr par-est
estf <- function(p){
  mu <- log(obs$x) + p["r"]*(1-obs$x/p["K"])
  -sum(dlnorm(obs$y, mu, p["s"]), log=TRUE)
}
o <- optim(par = c(r=1,K=mean(x),s=1), estf, method="L", lower=c(1e-3,1e-3,1e-3))
f_alt <- Ricker
p_alt <- c(o$par['r'], o$par['K'])
sigma_g_alt <- o$par['s']



## @knitr par-est-BH
estf <- function(p){
  mu <- obs$x * p["A"] / (1+ obs$x * p["B"])
  -sum(dlnorm(obs$y, log(mu), p["s"]), log=TRUE)
}
o <- optim(par = c(A=2, B=1/mean(x),s=1), estf, method="L", lower=c(1e-3,1e-3,1e-3))
f_alt <- BevHolt
p_alt <- c(A=as.numeric(o$par['A']), B=as.numeric(o$par['B']))
sigma_g_alt <- o$par['s']


## @knitr gp-priors
s2.p <- c(50,50) #inv gamma has mean b / (a - 1) (assuming a>1) and variance b ^ 2 / ((a - 2) * (a - 1) ^ 2) (assuming a>2)
tau2.p <- c(50,50)
d.p = c(10, 1/0.01, 10, 1/0.01)  ## sum of gammas, shape a & rate b. gamma has mean a / b and variance a / b ^ 2
nug.p = c(10, 1/0.01, 10, 1/0.01)
s2_prior <- function(x) dinvgamma(x, s2.p[1], s2.p[2])
tau2_prior <- function(x) dinvgamma(x, tau2.p[1], tau2.p[2])
d_prior <- function(x) dgamma(x, d.p[1], scale = d.p[2]) + dgamma(x, d.p[3], scale = d.p[4])
nug_prior <- function(x) dgamma(x, nug.p[1], scale = nug.p[2]) + dgamma(x, nug.p[3], scale = nug.p[4])
beta0_prior <- function(x, tau = tau2.p[2] / (tau2.p[1] - 1)) dnorm(x, 0, tau)
beta = c(0)


## @knitr gp-fit
gp <- bgp(X=obs$x, XX=x_grid, Z=obs$y, verb=0,
          meanfn="constant", bprior="b0", BTE=c(2000,16000,2),
          m0r1=FALSE, corr="exp", trace=TRUE, 
          beta = beta, s2.p = s2.p, d.p = d.p, nug.p = nug.p, tau2.p = tau2.p,
          s2.lam = "fixed", d.lam = "fixed", nug.lam = "fixed", tau2.lam = "fixed")


## @knitr gp-data
V <- gp$ZZ.ks2
Ef = gp$ZZ.km
tgp_dat <- data.frame(x   = gp$XX[[1]], 
                  y   = gp$ZZ.km, 
                 ymin = gp$ZZ.km - 1.96 * sqrt(gp$ZZ.ks2), 
                 ymax = gp$ZZ.km + 1.96 * sqrt(gp$ZZ.ks2),
                 ymin2 = gp$ZZ.mean - 1.96 * sqrt(gp$ZZ.vark), 
                 ymax2 = gp$ZZ.mean + 1.96 * sqrt(gp$ZZ.vark) )

# ZZ.ks2, ZZ.vark, ZZ.s2

## @knitr gp-plot
true <- sapply(x_grid, f, 0, p)
est <- sapply(x_grid, f_alt, 0, p_alt)
models <- data.frame(x=x_grid, GP=tgp_dat$y, Parametric=est, True=true)
models <- melt(models, id="x")
names(models) <- c("x", "method", "value")
# plot
ggplot(tgp_dat)  + 
  geom_ribbon(aes(x,y,ymin=ymin,ymax=ymax), fill="gray80") +
  geom_ribbon(aes(x,y,ymin=ymin2,ymax=ymax2), fill="gray60") +
  geom_line(data=models, aes(x, value, col=method), lwd=2, alpha=0.8) + 
  geom_point(data=obs, aes(x,y), alpha=0.8) + 
  xlab(expression(X[t])) + ylab(expression(X[t+1])) +
  scale_colour_manual(values=cbPalette)



## @knitr gp-posteriors
hyperparameters <- c("index", "s2", "tau2", "beta0", "nug", "d", "ldetK")
posteriors <- melt(gp$trace$XX[[1]][,hyperparameters], id="index")
priors <- list(s2 = s2_prior, tau2 = tau2_prior, beta0 = beta0_prior, nug = nug_prior, d = d_prior, ldetK = function(x) 0)
prior_curves <- ddply(posteriors, "variable", function(dd){
  grid <- seq(min(dd$value), max(dd$value), length = 100)
  data.frame(value = grid, density = priors[[dd$variable[1]]](grid))
})
ggplot(posteriors) + 
  #geom_density(aes(value), lwd=2) +
  geom_histogram(aes(x=value, y=..density..), lwd=2) +
  geom_line(data=prior_curves, aes(x=value, y=density), col="red", lwd=2) +
  facet_wrap(~ variable, scale="free")
ggplot(prior_curves) + 
  geom_line(aes(x=value, y=density), col="red", lwd=2) +
  facet_wrap(~ variable, scale="free")
#ggplot(subset(posteriors, variable=="nug")) + geom_histogram(aes(x=value, y = ..density..), lwd=2) + stat_function(fun = nug_prior, col="red", lwd=2)
#ggplot(subset(posteriors, variable=="s2")) + geom_histogram(aes(x=value, y = ..density..), lwd=2) + stat_function(fun = s2_prior, col="red", lwd=2)


## @knitr persistence-test
X <- numeric(length(x_grid))
X[38] = 1
h <- 0


## @knitr gp-F-sim 
F_ <- gp_F(h, Ef, V, x_grid)
xt1 <- X %*% F_
xt10 <- xt1
for(s in 1:OptTime)
  xt10 <- xt10 %*% F_
qplot(x_grid, xt10[1,]) + geom_point(aes(y=xt1[1,]), col="grey")



## @knitr par-F-sim
F_true <- par_F(h, f, p, x_grid, sigma_g)
yt1 <- X %*% F_true
yt10 <- yt1
for(s in 1:OptTime)
  yt10 <- yt10 %*% F_true
qplot(x_grid, yt10[1,]) + geom_point(aes(y=yt1[1,]), col="grey")


## @knitr F-sim-plot
transition <- melt(data.frame(x = x_grid, gp = xt1[1,], parametric = yt1[1,]), id="x")
ggplot(transition) + geom_point(aes(x,value, col=variable))



## @knitr est-F-sim
F_est <- par_F(h, f_alt, p_alt, x_grid, sigma_g)
zt1 <- X %*% F_est
zt10 <- zt1
for(s in 1:OptTime)
  zt10 <- zt10 %*% F_est
qplot(x_grid, zt10[1,]) + geom_point(aes(y=zt1[1,]), col="grey")




## @knitr gp-opt
matrices_gp <- gp_transition_matrix(Ef, .01*V, x_grid, h_grid)
opt_gp <- find_dp_optim(matrices_gp, x_grid, h_grid, OptTime, xT, profit, delta, reward=reward)

## @knitr true-opt
matrices_true <- f_transition_matrix(f, p, x_grid, h_grid, sigma_g)
opt_true <- find_dp_optim(matrices_true, x_grid, h_grid, OptTime, xT, profit, delta=delta, reward = reward)


## @knitr est-opt
matrices_estimated <- f_transition_matrix(f_alt, p_alt, x_grid, h_grid, sigma_g_alt)
opt_estimated <- find_dp_optim(matrices_estimated, x_grid, h_grid, OptTime, xT, profit, delta=delta, reward = reward)


## @knitr policy_plot
policies <- melt(data.frame(stock=x_grid, 
                            GP = x_grid[opt_gp$D[,1]], 
                            Parametric = x_grid[opt_estimated$D[,1]],
                            True = x_grid[opt_true$D[,1]]),
                  id="stock")
names(policies) <- c("stock", "method", "value")
policy_plot <- ggplot(policies, aes(stock, stock - value, color=method)) +
  geom_line(lwd=2, alpha=0.8) + 
  xlab("stock size") + ylab("escapement")  +
  scale_colour_manual(values=cbPalette)
policy_plot


## @knitr othernoise
z_g = function() rlnorm(1, 0, sigma_g)
z_m = function() 1+(2*runif(1, 0,  1)-1) * sigma_m

## @knitr stationary_policy_only
m <- sapply(1:OptTime, function(i) opt_gp$D[,1])
opt_gp$D <- m
mm <- sapply(1:OptTime, function(i) opt_true$D[,1])
opt_true$D <- mm
mmm <- sapply(1:OptTime, function(i) opt_estimated$D[,1])
opt_estimated$D <- mmm



## @knitr simulate
set.seed(1)
sim_gp <- lapply(1:100, function(i) ForwardSimulate(f, p, x_grid, h_grid, K, opt_gp$D, z_g, profit=profit))
set.seed(1)
sim_true <- lapply(1:100, function(i) ForwardSimulate(f, p, x_grid, h_grid, K, opt_true$D, z_g, profit=profit))
set.seed(1)
sim_est <- lapply(1:100, function(i) ForwardSimulate(f, p, x_grid, h_grid, K, opt_estimated$D, z_g, profit=profit))


## @knitr tidy
dat <- list(GP = sim_gp, Parametric = sim_est, True = sim_true)
dat <- melt(dat, id=names(dat[[1]][[1]]))
dt <- data.table(dat)
setnames(dt, c("L1", "L2"), c("method", "reps")) 


## @knitr sim-fish
ggplot(dt) + 
  geom_line(aes(time, fishstock, group=interaction(reps,method), color=method), alpha=.1) +
  scale_colour_manual(values=cbPalette, guide = guide_legend(override.aes = list(alpha = 1)))


## @knitr sim-harvest
ggplot(dt) +
  geom_line(aes(time, harvest, group=interaction(reps,method), color=method), alpha=.1) +
  scale_colour_manual(values=cbPalette, guide = guide_legend(override.aes = list(alpha = 1)))


## @knitr costs
profits <- dt[, sum(profit), by = c("reps", "method")]
means <- profits[, mean(V1), by = method]
sds <- profits[, sd(V1), by = method]
yield <- cbind(means, sd = sds$V1)
yield


