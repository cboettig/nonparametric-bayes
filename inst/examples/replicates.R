models <- c("Myers","Allen")

parameters <- list(Myers = list(
  c(r=1.5 + rnorm(1, 0, .1), theta=2.5 + rnorm(1, 0, .1), K=8 + rnorm(1, 0, .2)),
  c(r=1.5 + rnorm(1, 0, .1), theta=2.5 + rnorm(1, 0, .1), K=8 + rnorm(1, 0, .2)),
  c(r=1.5 + rnorm(1, 0, .1), theta=2.5 + rnorm(1, 0, .1), K=8 + rnorm(1, 0, .2))),
                   Allen = list(
  c(r=2 + rnorm(1, 0, .1), K=8 + rnorm(1, 0, .1), C=5 + rnorm(1, 0, .2)),
  c(r=2 + rnorm(1, 0, .1), K=8 + rnorm(1, 0, .1), C=5 + rnorm(1, 0, .2)),
  c(r=2 + rnorm(1, 0, .1), K=8 + rnorm(1, 0, .1), C=5 + rnorm(1, 0, .2)))
                   )
nuisance_pars <- c("sigma_g")
nuisance_values <- list(sigma_g = c(0.01, 0.05, 0.1))
replicates <- c(1111, 2222, 3333, 4444, 5555, 6666, 7777, 8888) # seeds


sensitivity <- function(model, parameters, nuisance, seed){
  
  if(model == "Myers")
    f <- Myers
  else if(model == "Allen")
    f <- RickerAllee
    
  sigma_g <- nuisance[["sigma_g"]]
  z_g <- function() rlnorm(1, 0, sigma_g)
  p <- parameters

  config = list(
    f = f,
    p = p,
    sigma_g = sigma.g,
    z_g = z_g,
    x_grid = seq(0, 15, length=50),
    h_grid = x_grid,
    profit = function(x,h) pmin(x, h),
    delta = 0.01,
    OptTime = 50,  # stationarity with unstable models is tricky thing
    reward = 0,
    xT = 0,
    Xo = 5.5,# observations start from
    x0 = 8, # simulation under policy starts from
    Tobs = 40,
    MaxT = 1000) # timeout for value iteration convergence

  attach(config)
  # replicate over random seed
  yields <- sapply(seed, 
                   function(seed_i){
                     set.seed(seed_i)

## Simulate data
                     x <- numeric(Tobs)
                     x[1] <- Xo
                     nz <- 1
                     for(t in 1:(Tobs-1))
                       x[t+1] = z_g() * f(x[t], h=0, p=p)
                     
                     X = c(rep(0,nz), pmax(rep(0,Tobs-1), x[1:(Tobs-1)]))
                     Y = c(rep(0,nz), x[2:Tobs])

## GP Stuff

                     ## @knitr gp-priors
                     s2.p <- c(5,5)  
                     d.p = c(10, 1/0.1)
                     
                     ## @knitr gp
                     gp <- gp_mcmc(X, y=Y, n=1e5, s2.p = s2.p, d.p = d.p)
                     gp_dat <- gp_predict(gp, x_grid, burnin=1e4, thin=300)
                     
                     matrices_gp <- gp_transition_matrix(gp_dat$Ef_posterior, gp_dat$Vf_posterior, x_grid, h_grid) 
                     opt_gp <- value_iteration(matrices_gp, x_grid, h_grid, MaxT, xT, profit, delta, reward)
                     
                     ## @knitr mle-opt
                     matrices_true <- f_transition_matrix(f, p, x_grid, h_grid, sigma_g)
                     opt_true <- value_iteration(matrices_true, x_grid, h_grid, OptTime=MaxT, xT, profit, delta=delta)
          
                     ## Assemble policies
                     OPT = data.frame(GP = opt_gp$D, True = opt_true$D)
                     policies <- melt(data.frame(stock=x_grid, sapply(OPT, function(x) x_grid[x])), id="stock")
                     names(policies) <- c("stock", "method", "value")
                     
                     
                     ## Simulations
                     sims <- lapply(OPT, function(D){
                       set.seed(1)
                       lapply(1:100, function(i) 
                         ForwardSimulate(f, p, x_grid, h_grid, x0, D, z_g, profit=profit, OptTime=OptTime)
                              )
                     })
                    
                     ## Assemble results
                     dat <- melt(sims, id=names(sims[[1]][[1]]))
                     sims_data <- data.table(dat)
                     setnames(sims_data, c("L1", "L2"), c("method", "reps")) 
                     # Legend in original ordering please, not alphabetical: 
                     sims_data$method = factor(sims_data$method, ordered=TRUE, levels=names(OPT))
                     
                     Profit <- sims_data[, sum(profit), by=c("reps", "method")]
                     tmp <- dcast(Profit, reps ~ method)
                     tmp <- tmp / tmp[,"True"]
                     tmp <- melt(tmp[2:dim(tmp)[2]])
                     tmp$value[tmp$variable == "GP"]  
                   })
 




  yields_dat <- melt(yields)
  names(yields_dat) <- c("replicate", "simulation", "value")
  
  # Make ids factors, not integers
  yields_dat$replicate <- as.factor(yields_dat$replicate)
  yields_dat$simulation <- as.factor(yields_dat$simulation)
  
  
  # definitions of id codes.  In this case the id number is it's own definition.  
  rep_ids <-  levels(yields_dat$replicate)
  names(rep_ids) <- rep_ids
  sim_ids <-  levels(yields_dat$simulation)
  names(sim_ids) <- sim_ids  
  
  
  
  dat <- data.frame(model = model, 
                    pars = as.list(parameters), 
                    replicate = yields_dat$replicate, 
                    sim = yields_dat$simulation, 
                    value = yields_dat$value, 
                    noise = sigma_g)
}
model <- "Allen"
allen1.01 <- sensitivity(model, 
                   parameters = parameters[[model]][[1]], 
                   nuisance = c(sigma_g = nuisance_values$sigma_g[1]), 
                   seed=c(1234, 2222, 3333))
model <- "Allen"
allen2.01 <- sensitivity(model, 
                   parameters = parameters[[model]][[2]], 
                   nuisance = c(sigma_g = nuisance_values$sigma_g[1]), 
                   seed=c(1234, 2222, 3333))
model <- "Allen"
allen1.05 <- sensitivity(model, 
                   parameters = parameters[[model]][[1]], 
                   nuisance = c(sigma_g = nuisance_values$sigma_g[2]), 
                   seed=c(1234, 2222, 3333))
model <- "Allen"
allen2.05 <- sensitivity(model, 
                   parameters = parameters[[model]][[2]], 
                   nuisance = c(sigma_g = nuisance_values$sigma_g[2]), 
                   seed=c(1234, 2222, 3333))

model <- "Myers"
Myers1.01 <- sensitivity(model, 
                   parameters = parameters[[model]][[1]], 
                   nuisance = c(sigma_g = nuisance_values$sigma_g[1]), 
                   seed=c(1234, 2222, 3333))
model <- "Myers"
Myers2.01 <- sensitivity(model, 
                   parameters = parameters[[model]][[2]], 
                   nuisance = c(sigma_g = nuisance_values$sigma_g[1]), 
                   seed=c(1234, 2222, 3333))
model <- "Myers"
Myers1.05 <- sensitivity(model, 
                   parameters = parameters[[model]][[1]], 
                   nuisance = c(sigma_g = nuisance_values$sigma_g[2]), 
                   seed=c(1234, 2222, 3333))
model <- "Myers"
Myers2.05 <- sensitivity(model, 
                   parameters = parameters[[model]][[2]], 
                   nuisance = c(sigma_g = nuisance_values$sigma_g[2]), 
                   seed=c(1234, 2222, 3333))
save(list=ls(), file="sensitivity.rda")


## Assemble into data.frame
allen_dat <- rbind(allen1.01, allen1.05, 
             allen2.01, allen2.05) 
m <- rbind(Myers1.01, Myers1.05, 
             Myers2.01, Myers2.05)
myers_dat <- m[c(1:2,4,3,5:8)]
names(myers_dat) <- names(allen_dat)
model_dat <- rbind(allen_dat, myers_dat)
dat <- model_dat
dat$pars.r <- factor(dat$pars.r, labels=c("A", "B", "C", "D"))
dat <- dat[c(1:2,5:6, 8, 7)]
dat$noise <- factor(dat$noise)
names(dat) <- c("model", "parameters", "replicate", "simulation", "noise", "value")

#ggplot(dat) + geom_histogram(aes(value)) + xlim(0,1.0) + theme_bw() + xlab("value as fraction of the optimal")
#ggplot(dat) + geom_histogram(aes(value, fill=noise)) + xlim(0,1.0) + theme_bw() + xlab("value as fraction of the optimal") 
#ggplot(dat) + geom_histogram(aes(value, fill=noise)) + xlim(0,1.0) + theme_bw() + xlab("value as fraction of the optimal") + facet_wrap(~model, ncol=1)





## Extract acutal parameter values corresponding to each parameter set
  p1 = levels(factor(model_dat$pars.r))
  p2 = levels(factor(model_dat$pars.K))
  p3 = levels(factor(model_dat$pars.C))

  A = c(r = p1[1], K = p2[1], theta = p3[1])
  B = c(r = p1[2], K = p2[2], theta = p3[2])
  C = c(r = p1[3], K = p2[3], C = p3[3])
  D = c(r = p1[4], K = p2[4], C = p3[4])

