Non-parametric approaches to optimal policy are more robust
===========================================================





Carl Boettiger, Steve Munch, Marc Mangel

# Abstract


# Introduction

The problem of structural uncertainty in managing ecological systems.  

* Most management recommendations from the ecological literature / management policies based on (motivated by) parametric models.  _Preference to pitch towards policy or theoretical literature?_

* Background discussion on the importance/success of parametric modeling; (e.g. @Levins1966, etc., up through @Geritz2011a) 

* Background on the concerns of structural uncertainty -- we don't have the right models.  (Also: measurement uncertainty, parameter uncertainty, unobserved states, boundary conditions, etc.)

* Hierarchical Bayesian approach has provided a natural way to address these from a statistical standpoint.  Successes and challenges from the parametric route: e.g. @Cressie2009.  

* Among pitfalls of these approaches: particularly difficult to apply in management context. Optimization-based (Decision-theoretic) approaches to ecological management need to be precisely parameterized about everything, even the uncertainty itself [@Polasky2011] (discussion of @Brozovic2011 SDP example in context of threshold system would perhaps be useful as well). 

* Management goals / decision-theoretic approaches need accurate prediction over relevant (short?) timescales more than accurate (but incomplete or noisy estimated) mechanisms.  


* Nonparametric (machine-learning?) approaches may offer the benefit of the hierarchical Bayesian approach without the practical and computational limitations of their parametric kin.  Non-parametric models are flexible enough to take maximum advantage of the data available, while being appropriately ambiguous about the dynamics of a system in regions of parameter space that have been poorly or never sampled.  

* Nonparametric approaches are beginning to appear more frequently in ecological and conservation literature (Species distribution models/maxent), including the Gaussian process based approach used here [@Munch2005].  However, such approaches have yet to be applied to the decision-theoretic framework that could guide management decisions, where we expect them to excel for several reasons (1) the ability to make accurate forecasts by more closely approximating the underlying process where data is available (2) remaining appropriately ambiguous where data is not available (3) remaining computationally simple enough to avoid some pitfalls common to hierarchical parametric approaches.  




# Approach and Methods







### Discussion of state equation




Concerns over the potential for tipping points in ecological dynamics [@Scheffer2001] highlight the dangers of uncertainty in ecological management and pose a substantial challenge to existing decision-theoretic approaches [@Brozovic2011].  To compare the performance of nonparametric and parametric approaches in an example that is easy to conceptualize, we will focus on a simple parametric model for a single species [derived from fist principles by @Allen2005a] as our underlying "reality".  

\begin{align}
X_{t+1} &= Z_t f(S_t) \\
S_t &= X_t - h_t \\
f(S_t) &= e^{r \left(1 - \frac{S_t}{K}\right)\left(S_t - C\right)}
\end{align}

As a low-dimensional system completely described by three parameters, this scenario should if anything be favorable to a parametric-based approach.  This model contains an Allee effect, or tipping point, below which the population is not self-sustaining and shrinks to zero [@Courchamp2008].





### Brief discussion on choice of model parameters, nuisance parameters

Where $Z_t$ is multiplicative noise function with mean 1, representing stochastic growth. We will consider log-normal noise with shape parameter $\sigma_g$.  We start with an example in which the parameters are $r =$ `2`, $K =$ `10`, $C =$ `5`, and  $\sigma_g =$ `0.05`.


### Discussion of training data

Both parametric and nonparametric approaches will require some training data on which to base their model of the process.  We generate the training data under the model described in Eq 1 for `35` time steps, under a known but not necessarily optimal sequence of harvest intensities, $h_t$.  For simplicity we imagine a fishery that started from zero harvest pressure and has been gradually increasing the harvest.  (Motivation, alternatives, stationarity, examples without a stable node (limit-cycle models), examples based on observations near a stable node alone, and why that isn't impossible).  







### Discussion of maximum likelihood estimated models




We estimate two parametric models from the data using a maximum likelihood approach.  The first model is structurally identical to the true model (Eq 1), differing only in that it's parameters are estimated from the observed data rather than given.  The alternative model is the Ricker model, which is structurally similar and commonly assumed 


(MLE models will assume the noise is log-normal, which it is in the simulation).  


Which estimates a Ricker model with $r =$ `1.8501`, $K =$ `9.8091`, and the Allen Allele model with $r =$ `2.8079`, $K =$ `11.8235` and $C =$ `7.2159`.  


### (Brief) Discussion of GP inference

* The use of Gaussian processes for inference in dynamical systems [introduced by @Kocijan2005]

* Gaussian processes in ecological literature [@Munch2005] 


* Our methodology (e.g. following @Munch2005).  








###  Discussion of the dynamic programming solution

(More thorough, but general-audience targeted.  Technical details and code provided in appendices).

Outside the ecological community, Gaussian processes have been introduced into optimization literature, but as an approximation to the value function rather than to underlying dynamics [@Deisenroth2009].  






### Discussion on how we compare performance of policies




* Replicate stochastic simulations 
* Sensitivity analysis (Figure 4).  

# Results

## Figure 1: 

_Shows the inferred Gaussian Process compared to the true and parametric models.  Refer to the appendix for details on the GP posteriors, etc._

![Graph of the inferred Gaussian process compared to the true process and maximum-likelihood estimated process.  Graph shows the expected value for the function $f$ under each model.  Two standard deviations from the estimated Gaussian process covariance with (light grey) and without (darker grey) measurement error are also shown.  The training data is also shown as black points.  (The GP is conditioned on 0,0, shown as a pseudo-data point). ](figure/gp_plot.pdf) 



## Figure 2: 

_The take-home message, showing that the GP is closest to the optimal strategy, while the parametric methods are less accurate.  Visualizing the policy may be more useful for the technical reader, the general audience may prefer Figure 3 showing all replicates of the population collapse under the parametric model and not under the GP._

![The steady-state optimal policy (infinite boundary) calculated under each model.  Policies are shown in terms of target escapement, $S_t$, as under models such as this a constant escapement policy is expected to be optimal [@Reed1979].](figure/policies_plot.pdf) 


## Figure 3: 

_Figure 3 is a less abstract and more visceral visualization of the take-home message, with the structurally inaccurate model leading universally to a collapse of the fishery and very few profits, while the Gaussian process performs nearly optimally.  The parametric approach even with the correct underlying structure does not perform optimally, choosing in this case to under-fish (may need to show harvest dynamics since that is not clear from the figure! Also isn't general, sometimes does optimally, sometimes over-fishes.  Perhaps need to show more examples.)  May need to show profits too?_ 

![Gaussian process inference outperforms parametric estimates. Shown are 100 replicate simulations of the stock dynamics (eq 1) under the policies derived from each of the estimated models, as well as the policy based on the exact underlying model.](figure/sim_plot.pdf) 


## Figure 4:

_Shows the sensitivity analysis.  A histogram of distribution of yield over stochastic realizations, showing that the qualitative results do not depend on the stochastic realization of the training data here, or on the parameters of the underlying model, though quantitative differences are visible._



# Discussion / Conclusion


* Non-parametric methods have received far too little attention in ecological modeling efforts that are aimed at improved conservation planning and decision making support.  

* Importance of non-parametric approaches in conservation planning / resource management / decision theory.  

* Decision-theoretic tools such as optimal control calculations rely on robust _forecasting_ more strongly than they rely on accurate _mechanistic_ relationships.  

* Adapting a non-parametric approach requires modification of existing methods for decision theory.  We have illustrated how this might be done in the context of stochastic dynamic programming, opening the door for substantial further research into how these applications might be improved.  

* Anticipate improved relative performance in higher dimensional examples

* Discuss constant escapement in model, in policies.  

* Limitations of this comparison: Are the maximum-likelihood solutions a straw man?

* Discussion of alternative related approaches: POMDP/MOMDP,  

## Future directions 

* Multiple species
* Online learning
* Multiple step-ahead predictions
* Explicitly accomidating additional uncertainties
* Improving inference of optimal policy from the GP





# Appendix / Supplementary Materials

## MCMC posterior distributions and convergence analysis

![Histogram of posterior distributions for the estimated Gaussian Process shown in Figure 1.  Prior distributions overlaid.](figure/posteriors.pdf) 

 
 @Gramacy2005
 
## Tables of nuisance parameters, sensitivity analysis

### List of hyper-parameters, prior distributions and their parameters

## Reproducible code, "Research Compendium"



