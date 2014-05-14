
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

opts_knit$set(fig.path="figure/small-allee/")

opts_chunk$set(fig.width=6, fig.height=6, echo=FALSE, dev="pdf", warning=FALSE, message=FALSE, echo=FALSE)
theme_set(theme_bw(base_size=12))
toggle = "hide" # results argument





load("small-allee.rda")

## ----figure_1, dependson=c("assemble-models", "par-fns"), fig.cap="Points show the training data of stock-size over time.  Curves show the posterior step-ahead predictions based on each of the estimated models. Observe that all models are fitting the data reasonably well.", fig.width=8, fig.height=6----
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


## ----figure_2, dependson=c("assemble-models"), fig.cap="Graph of the inferred Gaussian process compared to the true process and maximum-likelihood estimated process.  Graph shows the expected value for the function $f$ under each model.  Two standard deviations from the estimated Gaussian process covariance with (light grey) and without (darker grey) measurement error are also shown.  The training data is also shown as black points.  The GP is conditioned on (0,0), shown as a pseudo-data point.", fig.width=8, fig.height=6----
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


## ----figure_3, fig.width=8, fig.height=6, fig.cap="Beyond the tipping point (outside the range of the training data), the true dynamics fall outside the uncertainty bands of the wrong parametric models, but inside the uncertainty predicted by the GP. Points show the stock size simulated by the true model.  Overlay shows the range of states predicted by each model, based on the state observed in the previous time step. The Ricker model always (wrongly) predicts positive population growth, while the actual population shrinks in each step as the initial condition falls below the Allee threshold of the underlying model (Allen).  Note that the GP is both more pessimistic and more uncertain about the future state than the parametric models."----
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




## ----figure_S1, fig.cap = "Traces from the MCMC estimates of the GP model show reasonable mixing (no trend) and sampling rejection rate (no piecewise jumps)"----
gp_assessment_plots$traces_plot


## ----figure_S2, fig.cap="Posterior distributions from the MCMC estimate of the GP model. Prior curves shown in red."----
gp_assessment_plots$posteriors_plot


## ----figure_S3, fig.height=6, fig.cap="Traces from the MCMC estimates of the Ricker model show reasonable mixing (no trend) and sampling rejection rate (no piecewise jumps)", fig.width=6----
plot_ricker_traces


## ----figure_S4, fig.cap="Posteriors from the MCMC estimate of the Ricker model", fig.width=6, fig.height=4----
ggplot(ricker_posteriors, aes(value)) + 
  stat_density(geom="path", position="identity", alpha=0.7) +
  facet_wrap(~ variable, scale="free", ncol=2)


## ----figure_S5, fig.height=6, fig.cap="Traces from the MCMC estimates of the Myers model show reasonable mixing (no trend) and sampling rejection rate (no piecewise jumps)", fig.width=6----
plot_myers_traces


## ----figure_S6, fig.cap="Posterior distributions from the MCMC estimates of the Myers model", fig.width=6, fig.height=6----
ggplot(myers_posteriors, aes(value)) + 
  stat_density(geom="path", position="identity", alpha=0.7) +
  facet_wrap(~ variable, scale="free", ncol=2)


## ----figure_S7, fig.height=6, fig.cap="Traces from the MCMC estimates of the Allen model show reasonable mixing (no trend) and sampling rejection rate (no piecewise jumps)", fig.width=6----
plot_allen_traces


## ----figure_S8, fig.cap="Posteriors from the MCMC estimate of the Allen model", fig.width=6, fig.height=6----
ggplot(allen_posteriors, aes(value)) + 
  stat_density(geom="path", position="identity", alpha=0.7) +
  facet_wrap(~ variable, scale="free", ncol=2)


## ----include=FALSE-------------------------------------------------------
unlink("ricker_process.bugs")
unlink("allen_process.bugs")
unlink("myers_process.bugs")


save(list=ls(), file="small-allee.rda")
