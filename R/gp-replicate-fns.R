cbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")


#' Simulate observed data for training the GP
#' 
#' @param Xo the initial population size
#' @param z_g a function of no arguments providing the random multiplicative growth noise
#' @param f a function of (x, h, p) providing the growth in population x under harvest h given pars p
#' @param p the parameters of function f
#' @param Tobs the number of timesteps on which to observe the system
#' @param seed the random number seed for the simulation
#' @param nz the number of observations at 0 (a hack, should just condiition on 0,0)
#' @param harvest a sequence of harvest values that occur durning the simulation
#' @details varying harvest values allow the system to explore the state space, making for better training data.
#' @export
sim_obs <- function(Xo, z_g, f, p, Tobs = 35, seed = 1, nz = 10, 
                    harvest = sort(rep(seq(0, .5, length=7), 5)),
                    plotdata=FALSE){
  x <- numeric(Tobs)
  x[1] <- Xo
  set.seed(seed)
  for(t in 1:(Tobs-1))
    x[t+1] = z_g() * f(x[t], h=harvest[t], p=p)
  if(plotdata)
    plot(x)
  obs <- data.frame(x = c(rep(0,nz), 
                          pmax(rep(0,Tobs-1), x[1:(Tobs-1)] -
                                 harvest[1:Tobs-1])), 
                    y = c(rep(0,nz), 
                          x[2:Tobs]))
  #plot(obs$x, obs$y)
  obs
}


#' MLE estimate of parameters given by function f
#' @param obs the observed data, two columns giving x_t, x_t+1 respectively
#' @return a list with the f given, the MLE estimated parameters, and estimated noise level
#' @export 
par_est_allee <- function(obs, f, p, 
                          init = c(p[1]+1, p[2]-1, p[3]-2, 
                                   s = sigma_g + abs(rnorm(1,0,.1)))
                          ){
  estf <- function(p){ 
    mu <- f(obs$x,0,p)
    -sum(dlnorm(obs$y, log(mu), p[4]), log=TRUE)
  }
  par = init
  o <- optim(par, estf, method="L", lower=c(1e-3,1e-3,1e-3, 1e-3))
  f_alt <- f
  p_alt <- c(as.numeric(o$par[1]), as.numeric(o$par[2]), as.numeric(o$par[3]))
  sigma_g_alt <- as.numeric(o$par[4])
  list(f = f_alt, p = p_alt, sigma_g = sigma_g_alt)
}



#' MLE estimate of parameters given by function f
#' @param obs the observed data, two columns giving x_t, x_t+1 respectively
#' @return a list with the f given, the MLE estimated parameters, and estimated noise level
#' @export 
par_est <- function(obs,  init = c(r=1.5, K=mean(obs$x), s=1)){
  estf <- function(p){
    mu <- log(obs$x) + p["r"]*(1-obs$x/p["K"])
    llik <- -sum(dlnorm(obs$y, mu, p["s"]), log=TRUE)
    if(!is.numeric(llik) | is.nan(llik) | !(llik < Inf)){
      warning("possible error in llik")
      llik <- 1e30
    }
    llik
  }
  o <- optim(par = init, estf, method="L", lower=c(1e-3,1e-3,1e-3))
  f_alt <- Ricker
  p_alt <- c(o$par['r'], o$par['K'])
  sigma_g_alt <- o$par['s']
  list(f_alt = f_alt, p_alt = p_alt, sigma_g_alt = sigma_g_alt)
}


#' helper function to determine the optimal policies of each model
#' @import pdgControl
#' @export
optimal_policy <- function(gp, f, f_est, f_alt, p, p_est, p_alt, x_grid, h_grid, sigma_g, sigma_g_est, sigma_g_alt, delta, xT, profit, reward, OptTime){
  matrices_gp <- gp_transition_matrix(gp$ZZ.km, gp$ZZ.ks2, x_grid, h_grid)
  opt_gp <- find_dp_optim(matrices_gp, x_grid, h_grid, OptTime, xT, profit, delta, reward=reward)
  matrices_true <- f_transition_matrix(f, p, x_grid, h_grid, sigma_g)
  opt_true <- find_dp_optim(matrices_true, x_grid, h_grid, OptTime, xT, profit, delta=delta, reward = reward)
  matrices_estimated <- f_transition_matrix(f_est, p_est, x_grid, h_grid, sigma_g_est)
  opt_estimated <- find_dp_optim(matrices_estimated, x_grid, h_grid, OptTime, xT, profit, delta=delta, reward = reward)
  matrices_alt <- f_transition_matrix(f_alt, p_alt, x_grid, h_grid, sigma_g_alt)
  opt_alt <- find_dp_optim(matrices_alt, x_grid, h_grid, OptTime, xT, profit, delta=delta, reward = reward)
  
  gp_D <- sapply(1:OptTime, function(i) opt_gp$D[,1])
  true_D <- sapply(1:OptTime, function(i) opt_true$D[,1])
  est_D <- sapply(1:OptTime, function(i) opt_estimated$D[,1])
  alt_D <- sapply(1:OptTime, function(i) opt_alt$D[,1])
  
  list(gp_D = gp_D, true_D = true_D, est_D = est_D, alt_D = alt_D)
}

#' helper function to simulate the optimal policies of each model
#' @import pdgControl data.table reshape2
#' @export
simulate_opt <- function(OPT, f, p, x_grid, h_grid, x0, z_g, profit){
  set.seed(1)
  sim_gp <- lapply(1:100, function(i) ForwardSimulate(f, p, x_grid, h_grid, x0, OPT$gp_D, z_g, profit=profit))
  set.seed(1)
  sim_true <- lapply(1:100, function(i) ForwardSimulate(f, p, x_grid, h_grid, x0, OPT$true_D, z_g, profit=profit))
  set.seed(1)
  sim_est <- lapply(1:100, function(i) ForwardSimulate(f, p, x_grid, h_grid, x0, OPT$est_D, z_g, profit=profit))
  set.seed(1)
  sim_alt <- lapply(1:100, function(i) ForwardSimulate(f, p, x_grid, h_grid, x0, OPT$alt_D, z_g, profit=profit))
  
  dat <- list(GP = sim_gp, Parametric = sim_est, True = sim_true, Structural = sim_alt)
  dat <- melt(dat, id=names(dat[[1]][[1]]))
  dt <- data.table(dat)
  setnames(dt, c("L1", "L2"), c("method", "reps")) 
  setkey(dt, method) # chance the ordering
  dt[c("GP", "Parametric", "True", "Structural")]
}


#' helper function to draw the plots of fish stock and harvest dynamics
#' @param dt the data.table from simulate_opt
#' @import ggplot2
#' @return plots of stock and harvest dynamics over time and replicates
#' @export
sim_plots <- function(dt, seed=1, harvest=FALSE){
  fish_plot <- ggplot(dt) + 
    geom_line(aes(time, fishstock, group=interaction(reps,method), color=method), alpha=.1) +
    scale_colour_manual(values=cbPalette, guide = guide_legend(override.aes = list(alpha = 1)))+
    ggtitle(paste("seed", seed))
  print(fish_plot)
  harvest_plot <- ggplot(dt) +
    geom_line(aes(time, harvest, group=interaction(reps,method), color=method), alpha=.1) +
    scale_colour_manual(values=cbPalette, guide = guide_legend(override.aes = list(alpha = 1))) +
    ggtitle(paste("seed", seed))
  if(harvest)
    print(harvest_plot)
}

#' function to provide summary stats or the profits
#' @param dt the data table from simulate_opt
#' @return a table of mean and variance in profit over the replicates, by method
#' @import data.table
#' @export
profits_stats <- function(dt){
  profits <- dt[, sum(profit), by = c("reps", "method")]
  means <- profits[, mean(V1), by = method]
  sds <- profits[, sd(V1), by = method]
  yield <- cbind(means, sd = sds$V1)
  yield
}


#' plot the gaussian process, true model, and fitted parametric model(s)
#' @import ggplot2 reshape2 
#' @export
gp_plot <- function(gp, f, p, f_est, p_est, f_alt, p_alt, x_grid, obs, seed){
  tgp_dat <- 
    data.frame(  x = gp$XX[[1]], 
                 y = gp$ZZ.km, 
                 ymin = gp$ZZ.km - 1.96 * sqrt(gp$ZZ.ks2), 
                 ymax = gp$ZZ.km + 1.96 * sqrt(gp$ZZ.ks2),
                 ymin2 = gp$ZZ.mean - 1.96 * sqrt(gp$ZZ.vark), 
                 ymax2 = gp$ZZ.mean + 1.96 * sqrt(gp$ZZ.vark))
  true <- sapply(x_grid, f, 0, p)
  alt <- sapply(x_grid, f_alt, 0, p_alt)
  est <- sapply(x_grid, f_est, 0, p_est)
  models <- data.frame(x=x_grid, GP=tgp_dat$y, Parametric=est, True=true, Structural=alt)
  models <- melt(models, id="x")
  names(models) <- c("x", "method", "value")
  plot_gp <- ggplot(tgp_dat) + geom_ribbon(aes(x,y,ymin=ymin,ymax=ymax), fill="gray80") +
    geom_ribbon(aes(x,y,ymin=ymin2,ymax=ymax2), fill="gray60") +
    geom_line(data=models, aes(x, value, col=method), lwd=2, alpha=0.8) + 
    geom_point(data=obs, aes(x,y), alpha=0.8) + 
    xlab(expression(X[t])) + ylab(expression(X[t+1])) +
    scale_colour_manual(values=cbPalette) +
    ggtitle(paste("seed", seed))
  print(plot_gp)
}


#' plot the gp posteriors
#' @param gp a gp fit
#' @param priors prior distribution functions so that they can be shown on plot as well
#' @return the posteriors plot as histogram, with priors overlaid
#' @export
posteriors_plot <- function(gp, priors){
  hyperparameters <- c("index", "s2", "tau2", "beta0", "nug", "d", "ldetK")
  posteriors <- melt(gp$trace$XX[[1]][,hyperparameters], id="index")
  prior_curves <- ddply(posteriors, "variable", function(dd){
    grid <- seq(min(dd$value), max(dd$value), length = 100)
    data.frame(value = grid, density = priors[[dd$variable[1]]](grid))
  })
  plot_posteriors <- ggplot(posteriors) + 
    #geom_density(aes(value), lwd=2) +
    geom_histogram(aes(x=value, y=..density..), alpha=0.7) +
    geom_line(data=prior_curves, aes(x=value, y=density), col="red", lwd=2) +
    facet_wrap(~ variable, scale="free")
  print(plot_posteriors)
}




#' policy plot
#' @param x_grid stock grid
#' @param gpD gaussian process policy
#' @param estD an estimated policy (e.g. parametric uncertainty only)
#' @param trueD the true optimal policy
#' @param altD an alternative (structurally incorrect) policy
#' @return the policy plot
#' @export
plot_policies <- function(x_grid, gpD, estD, trueD, altD){
  policies <- melt(data.frame(stock=x_grid, 
                              GP = x_grid[gpD], 
                              Parametric = x_grid[estD],
                              True = x_grid[trueD],
                              Structural = x_grid[altD]),
                   id="stock")
  names(policies) <- c("stock", "method", "value")
  policy_plot <- ggplot(policies, aes(stock, stock - value, color=method)) +
    geom_line(lwd=2, alpha=0.8) + 
    xlab("stock size") + ylab("escapement")  +
    scale_colour_manual(values=cbPalette)
  print(policy_plot)
}
