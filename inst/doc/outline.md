Abstract
========

Introduction
============

Briefly set the stage: (see related papers)

- Decision-theory / optimal control framework
- Fisheries context 

@Fischer 
@Sethi2005, @Peretti, @Munch, 


(Classical division)

- Dynamic uncertainty (stochasticity)
- Parametric uncertainty
- Structural uncertainty




## Two central challenges: 

- We don't know the model
- We don't have data from where we need it most

## Why we don't know the model

- Complex dynamics
- model choice and model averaging approaches

## Why we don't have data where we need it

- Concerns of tipping points 
- Danger of learning 


- Note on "not magic": honest uncertainty + SDP
- Note on comparing models (via value function rather than by "fit")
- Note on terminology: "Non-parametric" 

## Background on the Gaussian Process

- Definition
- Previous application
- Why it is particularly suited to these two problems
- (Why this is a novel application thereof)

Approach and Methods
====================

### Summary of approach

### Statement of the optimal control problem

- Underlying model 
- Available data 
- Value function



For simplicity we assume profit is simply linear in the realized harvest (only
enforcing the restriction that harvest can not exceed available stock)


### Parametric models

- Statement of the models

We consider three candidate parametric models of the stock-recruitment
dynamics: The Ricker model, the Allen model [Allen 2005](), the Myers
model. The familiar Ricker model involves two parameters, corresponding
to a growth rate and a carrying capacity, and cannot support alternative
stable state dynamics (though as growth rate increases it exhibits a
periodic attractor that proceeds through period-doubling into chaos. We
will generally focus on dynamics below the chaotic threshold for
the purposes of this analysis.) The Allen model resembles the Ricker
dynamics with an added Allee effect parameter [Courchamp](), below
which the population cannot persist.  The Myers model also has three
parameters and contains an Allee threshold, but has compensatory rather
than over-compensatory density dependence (resembling a Beverton-Holt
curve rather than a Ricker curve at high densities.)

We assume multiplicative log-normal noise perturbs the growth predicted 
by the each of the deterministic model skeletons described above. This 
introduces one additional parameter $\sigma$ that must be estimated by each
model. 

<!-- equations just in appendix? -->

As we simulate training data from the Allen model (ref section), we will
refer to this as the structurally correct model.  The Ricker model is
thus a reasonable approximation of these dynamics far from the Allee
threshold (but lacks threshold dynamics), while the Myers model shares
the essential feature of a threshold but differs in the structure. Thus
we have three potential parametric models of the stock dynamics.

- Bayesian inference of parametric models

We infer posterior distributions for the parameters of each model
in a Bayesian context using Gibbs sampling (implemented in R [R]()
using jags, [jags]()).  We choose uninformative uniform priors for all
parameters (See Appendix, Figures S1-S3, and Table S1, and the R code
provided). One-step-ahead predictions of these model fits are shown in
Figure 1.

- SDP via parametric models

An optimal policy function is then inferred through stochastic dynamic
programming for each model given the posterior distributions of the
parameter estimates.  This policy maximizes the expectation of the value
function integrated over the parameter uncertainty. (code implementing
this algorithm provided in the Appendix).


### The Gaussian Process model

- Statement of model

... more on GP ... [Munch 2005]()

We also estimate a simple Gaussian Process defined by
a radial basis function kernel of two parameters: $\ell$, which gives
the characteristic length-scale over which correlation between two 
points (e.g. any two points $X_t, X_{t+1}$, and $X_{t+\tau}, X_{t+1+\tau}$)
in state-space decays, and $\sigma$, which gives the scale of the 
process noise by which observations $Y_{t+1}$ may differ from their
predicted values $X_{t+1}$ given an observation of the previous state,
$X_t$. 

- Inference of the model

Also unlike parametric models, this posterior distribution is still
conditional on the training data. As such, the uncertainty near the
observed data.

We use a Metropolis-Hastings Markov Chain Monte Carlo to infer posterior
distributions of the two parameters of the GP (Figure S4, code in
appendix), under weakly informative Gaussian priors (see parameters in
table S5). As the posterior distributions differ substantially from the
priors (Figure S4), we can be assured that most of the information in
the posterior comes from the data rather than the prior belief.


- SDP via the model

Though we are unaware of prior application of this type, it is reasonably
straight-forward to adapt the Gaussian Process for Stochastic Dynamic
Programming.  Recall that unlike the parametric models the Gaussian
process with fixed parameters already predicts a distribution of
curves rather than a single curve. We must first integrate over 
this distribution of curves given a sampling of parameter values drawn
from the posterior distribution of the two GP parameters, before
integrating over the posterior of those parameters themselves.



Results
=======

### Figure 1: Fitted Models

- All models fit the data quite well
- Information criteria would pick the simple, incorrect model.

(State-space figure??)

### Figure 2: Inferred Policies

- Inferred policies differ substantially among models
- The structurally correct model and the GP are close to the true model
- alternatives are not close

### Figure 3: Simulated results

- 

### Figure 4: Robustness

- Results hold across range of parameters

(eek, distinguish better between result and discussion?)

Discussion 
==========

- All models are "good fits" to the originally observed data. 
- (Simple model choice immediately leads us astray)


Future directions
-----------------

### Higher dimensions 

### Online learning

Acknowledgments
===============

Appendix
========

## Model definitions and estimation

Equation S1: Ricker model.

$$X_{t+1} = Z_t X_t e^{r \left(1 - \frac{S_t}{K} \right) } $$

Figure S1: Ricker model: prior and posterior distributions for parameter estimates.

```{r, echo=FALSE}
plot_ricker_traces
```

```{r, echo=FALSE}
plot_ricker_posteriors
```
Table S1: Parameterization of the priors

```{r echo=FALSE, results="asis"}
xtable::xtable(data.frame(
  parameter = c("r0", "K", "sigma"),
  lower_bound = c(r0_prior_p[1], K_prior_p[1], stdQ_prior_p[1]),
  upper_bound = c(r0_prior_p[2], K_prior_p[2], stdQ_prior_p[2])))
```


$$ X_{t+1} = Z_t \frac{r S_t^{\theta}}{1 - \frac{S_t^\theta}{K}} $$

Eq S2: Myers model 
Figure S2: Myers model: Traces, prior and posterior distributions for parameter estimates.

```{r, echo=FALSE}
plot_myers_traces
```

```{r, echo=FALSE}
plot_myers_posteriors
```

Table S2: Parameterization of the priors
```{r echo=FALSE, results="asis"}
xtable::xtable(data.frame(parameter = c("r0", "K", "theta", "sigma"),
           lower_bound = c(r0_prior_p[1], K_prior_p[1], theta_prior_p[1], stdQ_prior_p[1]),
           upper_bound = c(r0_prior_p[2], K_prior_p[2], theta_prior_p[2], stdQ_prior_p[2])))
```


Eq S3: Allen model 

$$f(S_t) = S_t e^{r \left(1 - \frac{S_t}{K}\right)\left(S_t - C\right)} $$

Figure S3: Allen model: prior and posterior distributions for parameter estimates.

```{r, echo=FALSE}
plot_allen_traces
```

```{r, echo=FALSE}
plot_allen_posteriors
```

Table S3: Parameterization of the priors

```{r echo=FALSE, results="asis"}
xtable::xtable(data.frame(parameter = c("r0", "K", "theta", "sigma"),
           lower_bound = c(r0_prior_p[1], K_prior_p[1], theta_prior_p[1], stdQ_prior_p[1]),
           upper_bound = c(r0_prior_p[2], K_prior_p[2], theta_prior_p[2], stdQ_prior_p[2])))
```



Eq S4: GP model 
Figure S4: GP model: prior and posterior distributions for parameter estimates.

```{r echo=FALSE}
gp_assessment_plots
```


Table S4: Parameterization of the priors

## Optimal Control Problem

We seek the harvest policy $h(x)$ that maximizes:

$$ \max_{h(t)} \sum_{t \in 0}^{\infty} e^{-\rho t} \Pi_t(X_t, h) \delta^t  $

subject to the profit function $\Pi(X_t,h)$, discount rate $\delta$, and the state
equation

$$X_{t+1} = Z_t f(S_t)  $$
$$S_t = X_t - h_t $$


Where $Z_t$ is multiplicative noise function with mean 1, representing
stochastic growth. We will consider log-normal noise with shape parameter
$\sigma_g$. 


Form this we can write down the Bellman recursion as: 

$$V_t(x_t) = \max_h \mathbf{E} \left(\Pi(h_t, x_t) + \delta V_{t+1}( Z_{t+1} f(x_t - h_t)) \right)$$

For simplicity we assume profit is simply linear in the realized harvest (only
enforcing the restriction that harvest can not exceed available stock), $\Pi(h,x) = min(h,x)$. 


### Pseudocode for the Bellman iteration

```r
 V1 <- sapply(1:length(h_grid), function(h){
      delta * F[[h]] %*% V +  profit(x_grid, h_grid[h]) 
    })
    # find havest, h that gives the maximum value
    out <- sapply(1:gridsize, function(j){
      value <- max(V1[j,], na.rm = T) # each col is a diff h, max over these
      index <- which.max(V1[j,])  # store index so we can recover h's 
      c(value, index) # returns both profit value & index of optimal h.  
    })
    # Sets V[t+1] = max_h V[t] at each possible state value, x
    V <- out[1,]                        # The new value-to-go
    D[,OptTime-time+1] <- out[2,]       # The index positions
```





### Training data

Eacho of our models $f(S_t)$ must be estimated from training data, which
we simulate from the Allen model with parameters $r = 2$,
$K=8$, $C=5$ and $\sigma_g = 0.1$ for $T = 40$ timesteps, starting 
at initial condition $X_0 = $.

<!-- ` r p[1]`, $K =$ ` r p[2]`, $C =$ ` r p[3]`, and  $\sigma_g =$ ` r sigma_g`. -->



