<!-- OLD TEXT -->

Abstract
=======================================================================


Decision-theoretic methods often rely on simple parametric models of
ecological dynamics to compare the value of a potential sequence of
actions. Unfortunately, such simple models rarely capture the complexity
or uncertainty found in most real ecosystems. 

Further, the data on which a model has been parameterized frequently 
fails to cover the possible state-space over which management decisions 
must operate.  Consequently a model do well in the region of state-space
in which it was estimated, but give erroneous confidence to predictions
outside of that region.  

This problem is keenly felt in any system where a potential threshold
or tipping point is a concern.  Such a tipping point, if it exists 
at all, will lay outside the observed range of the observed data. 


We demonstrate how nonparametric Bayesian models can provide robust,
solutions to decision making under uncertainty without knowing the 
structural form of the true model.  

While methods that account for _parametric_ uncertainty can be very 
successful with the right model,
structural uncertainty of not knowing what model best approximates the 
dynamics poses considerably greater difficulty.  



Introduction
=======================================================================

#### Opening 

<!-- More on complex dynamics and not having the correct model 


-->

<!-- More on the lack of data throughout the relevant state-space
     and how any concern about potential tipping points indicates 
     that the data exhibits this bias / problem.  
     -->




#### Models for decision-making under uncertainty 

Decision-theoretic or optimal control tools require a model that can assign probabilities of
future states (e.g. stock size of a fishery) given the current state and a
proposed action (e.g. fishing harvest or effort).  
Management frequently faces a sequential decision-making problem -- after 
selecting an action, the decision-maker may receive new information about 
the current state and must again choose an appropriate action -- such as setting
the harvest limits each year based on stock assessments the year prior.  

The decision maker typically seeks to determining the course of actions (also referred to as the policy)
that maximizes the expected value of some objective function such as
net present value derived from the resource over time.  
Though much can be said on how to choose this value function appropriately 
(e.g. see [@Halpern2013](http://doi.org/10.1073/pnas.1217689110 
"Achieving the triple bottom line in the face of inherent trade-offs 
among social equity, economic return, and conservation.")) we will
assume this is given.  (Nor is this approach necessarily constrained to 
maximizing the expectated value of such a function - the decision-theoretic
framework can be adapted to alternatives such as minimizing the maximum 
cost or damage that might be incurred; see @Polasky2011).

In representing future states with probabilities and maximizing expectations,
this approach provides a natural framework for handling uncertainty. 


The value function typically depends on the action or policy taken, as well as
the state of the system, in each interval of time. The state of the system,
in turn, is usually described by a dynamical model.


[@Williams2001; @Athanassoglou2012]. 


<!-- Transition and map: The weakness of parametric models -->
<!-- This section is a bit weak and vague.  Revisit to make a real map in light 
of outline structure below. -->

While simple mechanistic models can nevertheless provide important insights
into long-term outcomes, such approaches are not well-suited for use in
forecasting outcomes of potential management options.  Non-parametric
approaches offer a more flexible alternative that can both more accurately
reflect the data available while also representing greater uncertainty
in areas (of state-space) where data is lacking.

We demonstrate how
a Gaussian Process model of stock recruitment can lead to nearly optimal
management through stochastic dynamic programming, comperable to knowing
the correct structural equation for the underlying simulation.  Meanwhile,
parametric models that do not match the underlying dynamics can perform 
very poorly, even though they fit the data as well as the true model.  
Ecological research and management strategy should pay closer attention
to the opportunities and challenges nonparametric modeling can offer.



<!-- 
### Quantitative vs Qualitative Decisions

In this paper, we consider those ecological management problems in which
a mathematical (or computational) model is used to quantitatively inform
decision-making by allowing a manager to compare to the expected consequences
of potential management actions (or policies).  We distinguish this from the
solely qualitative use of a model, in which models are used to represent and
compare hypotheses of different mechanisms that would lead to qualitatively 
different actions. 

In addition to facing the potentially grave consequences of such transitions,
if they do exist in a system of interest, this generally means that we lack
data in the region near and beyond a possible tipping point.  As a consequence,
our models are least accurate just where we need them to be most accurate.  
Unfortunately, parametric modeling approaches are not well-suited to this 
scenario, in which good data is avialble over only part of the relevant 
state-space.  Necessary assumptions about model structure in any parametric
approach can lead such models to perform very poorly, even when accounting
for parameter uncertainty. 

-->


Approach and Methods
====================

## The optimal control problem in fisheries management

We focus on the problem in which a manager must set
the harvest level for a marine fishery each year to maximize the net
present value of the resource, given an estimated stock size from the
year before. 


<!-- ugh, re-word this -->

To permit comparisons against a theoretical optimum we will consider
data on the stock dynamics simulated from a simple parametric model
in which recruitment of the fish stock $X_{t+1}$ in the following year
is a stochastic process governed by a function $f$ of the current 
stock $X_t$, selected harvest policy $h_t$, and noise process $Z$,

$$X_{t+1} = Z_t f(X_t, h_t) $$


Given parameters for the function $f$ and probability distribution $Z$,
along with a given economic model determining the  price/profit $\Pi(X_t,
h_t)$ realized in a given year given a choice of harvest $h_t$ and
observed stock $X_t$.  This problem can be solved exactly for discretized
values of stock $X$ and policy $h$ using stochastic dynamic programming
(SDP) [@Mangel1988]. Problems of this sort underpin much marine fisheries
management today.

A crux of this approach is correctly specifying the functional form of $f$,
along with its parameters.  The standard approach uses one of a 
handful of common parametric models representing the stock-recruitment
relationship, usually after estimating the model parameters from any 
available existing data. Uncertainty in the parameter estimates can 
be estimated and integrated over to determine the optimal policy under
under uncertainty [@Mangel1988; @Schapaugh2013]. Uncertainty in the model
structure itself can only be addressed in this approach by hypothesizing 
alternative model structures, and then performing some model choice or
model averaging  [@Williams2001; @Athanassoglou2012]. 


## Underlying Model

To illustrate the value of the non-parametric Bayesian approach to management,
we focus on example of a system containing such a tipping point whose dynamics
can still be described by a simple, one-dimensional parametric model.  
We will focus on a simple parametric model for a single species [derived
from fist principles by @Allen2005a] as our underlying "reality".

$$X_{t+1} = Z_t f(S_t)  $$
$$S_t = X_t - h_t $$
$$f(S_t) = S_t e^{r \left(1 - \frac{S_t}{K}\right)\left(S_t - C\right)} $$

Where $Z_t$ is multiplicative noise function with mean 1, representing
stochastic growth. We will consider log-normal noise with shape parameter
$\sigma_g$.  We start with an example in which the parameters are $r =2$,
$K=8$, $C=5$ and $\sigma_g = 0.1$.


<!-- ` r p[1]`, $K =$ ` r p[2]`, $C =$ ` r p[3]`, and  $\sigma_g =$ ` r sigma_g`. -->


As a low-dimensional system completely described by three parameters, this
scenario should if anything be favorable to a parametric-based approach.
This model contains an Allee effect, or tipping point, below which the
population is not self-sustaining and shrinks to zero [@Courchamp2008].


#### Simulated training data

We generate initial observational data under the model described in Eq 1 
for $T_{\textrm{obs}}=40$ time steps, under a given arbitrary sequence of harvest
intensities, $h_t$. 
We consider the case in which most of the data comes from a limited region
of state space (e.g. near a stable equilibrium), leaving us without
observations of the population dynamics at very low levels which would
be useful in discrimating between recruitment curves [@] or demonstrating
the existence of a tipping point [@Scheffer2001].  
Using data simulated from a specified model rather than empirical data
permits the comparison against the true underlying dynamics, setting 
a bar for the optimal performance possible.  



### Parametric Models

We consider three candidate parametric models for the stock-recruitment function, 
which we refer to by the first authors of the publications in which they were
first proposed.  

We generate the data with a four-parameter model that contains
a tipping point, as discussed above (equation 1),
(an Allee effect, see [@Allen, @Courchamp]) below which the stock
decreases to zero, 

$$ X_{t+1} = Z_t S_t e^{r \left(1 - \frac{S_t}{K}\right)\left(\frac{S_t - \theta}{K}\right)} $$

$$ S_t = X_t - h_t $$

The parameter $C$ reflects the location of the tipping point, $K$ the carrying capacity of
the stock, and $r$ the base recruitment rate.  $S_t$ represents the stock size after
a harvest $h_t$ has been implemented.  $Z_t$ represents a log-normal random variable of
log-mean zero and log-standard deviation parameter $\sigma$.  


We consider two alternative candidate models: the Ricker [@Ricker] stock-recruitment curve,

$$X_{t+1} = Z_t X_t e^{r \left(1 - \frac{S_t}{K} \right) } $$

and an alternative four-parameter model adapted from @Myers, 

$$ X_{t+1} = Z_t \frac{r S_t^{\theta}}{1 - \frac{S_t^\theta}{K}} $$

which contains a tipping point for $\theta > 2$ and becomes a Beverton-Holt model 
at $\theta = 1$.  

<!-- Note that while similarly-named parameters contain similar connotations, the 
estimate and precise intepretation will be different in each model. -->


### Bayesian Inference of Parametric models 

Given the sample data, we infer posterior distributions for each of the three
models listed above using a Markov Chain Monte Carlo Gibbs Sampler (jags, see 
appendix for implementation details and code) given uniform priors.  We run six
chains for $10^6$ steps each and then assess convergence by Gelman-Rubin criterion
and inspection of the traces, see appendix.  

<!-- By transforming variables we can integrate out the the $\sigma$ parameter
and the constant term in the log-mean, see appendix.  The resulting one or
two parameters can then be evaluated over a vector or grid of possible values,
and the result transformed back to determine the posterior density of the parameters,
This approach avoids any concerns of convergence in the MCMC analysis and can thus
be used to verify that the performance of the parametric models can not be explained
by convergence issues.  However, this approach does not generalize easily to richer models,
such as the introduction of measurement noise or additional parameters, in which 
MCMC based approaches are the obvious choice.  -->




### The Non-parametric Bayesian alternative for stock-recruitment curves




### SDP via GP 

Once the posterior Gaussian process (GP) has been estimated [e.g. see
@Munch2005], it is necessary to adapt it in place of the parametric
equation for the stochastic dynamic programming (SDP) solution [see
@Mangel1988 for a detailed description of parametric SDP methods] to the
optimal policy. The essence of the idea is straight forward -- we will use
the estimated GP in place of the parametric growth function to determine
the stochastic transition matrix on which the SDP calculations are based.
The SDP is solved in a discretized state space -- both the continuously 
valued population densities $X$ and harvest quotas $h$ are first mapped to
a bounded, discrete grid.  (For simplicity we will consider a uniform grid,
though for either parametric or GP-based SDP it is often advantageous to 
use a non-uniform discretization such as a basis function representation,
e.g. see [@Deisenroth2009]). 

The SDP approach then computes a transition matrix, $\mathbf{F}$.  We 
demonstrate that calculation is just as straight forward based on the GP
as it is in the classical context using the parametric model.  The 
${i,j}$ of the transition matrix $F$ entry gives the probability of transitioning into 
state $x_i$ given that the system is in state $x_j$ in the previous 
time-step.  To generate the transition
matrix based on the posterior GP, we need only the expected values
at each grid point and the corresponding variances (the diagonal of
the covariance matrix), as shown in Figure 1.  Given the mean of the 
GP posterior at each grid-point as the vector $E$ and variance at that
point as vector $V$, the probability of transitioning from state $x_i$ to state $x_j$ is

$$\mathcal{N}\left(x_j | \mu = E_i, \sigma = \sqrt{V_i}\right)$$

where $\mathcal{N}$ is the Normal density at $x_j$ with mean $\mu$ and
variance $\sigma^2$.  Strictly speaking, the transition probability should
be calculated by integrating the normal density over the bin of width
$\Delta$ centered at $x_j$.  For a sufficiently fine grid that $f(x_j)
\approx f(x_j + \Delta)$, it is sufficient to calculate the density at
$x_j$ and then row-normalize the transition matrix. The process
can then be repeated for each possible discrete value of our control 
variable, (harvest $h$).  


**Pseudocode for the determining the transition matrix from the GP**

```r
for(h in h_grid)
  F_h = for(x_j in grid)
          for(i in 1:N) 
            dnorm(x_j, mu[i]-h, V[i])
```


Using the discrete transition matrix we may write down the Bellman
recursion defining the stochastic dynamic programming iteration:

\begin{equation}
V_t(x_t) = \max_h \mathbf{E} \left( h_t + \delta V_{t+1}( Z_{t+1} f(x_t - h_t)) \right)
\end{equation}

where $V(x_t)$ is the value of being at state $x$ at time $t$, $h$
is control (harvest level) chosen. Numerically, the maximization is
accomplished as follows. Consider the set of possible control values to
be the discrete values corresponding the grid of
stock sizes.  Then for each $h_t$ there is a corresponding transition
matrix $\mathbf{F}_h$ determined as described above but with mean 
$\mu = x_j - h_t$. Let $\vec{V_t}$ be the vector whose $i$th element corresponds
to the value of having stock $x_i$ at time $t$.  Then let $\Pi_h$ be
the vector whose $i$th element indicates the profit from harvesting
at intensity $h_t$ given a population $x_i$ (e.g. $\max(x_i, h_t)$
since one cannot harvest more fish then the current population size).
Then the Bellman recursion can be given in matrix form as

$$V_{t} = \max_h \left( \Pi_{h_{t}} + \delta \mathbf{F}_h V_{t+1} \right)$$

where the sum is element by element and the expectation is computed by the matrix multiplication $\mathbf{F} V_{t+1}$.  

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

This completes the algorithm adapting the GP to the sequential decision-making 
problem through SDP, which has not been previously demonstrated.  
We further provide an R package implementation as described in the supplemental materials.  


### Estimating parametric models 

We estimate posterior distributions for two parametric models: one using the structurally correct
model as given in Eq (1), which we refer to as the "Parametric Bayes" model,
and another using the familiar Ricker model, using a Gibbs sampler
as described (with source code) in the appendix).  In addition we 
estimate the parameters of the structurally correct model by maximum
likelihood.  
  

Results
=======

Discussion 
==========


#### Big picture: Linking GP to SDP  

_rambling_

Non-parametric Bayesian methods have received far too little attention
in ecological modeling efforts that are aimed at improved conservation
planning and decision making support.  Such approaches may be particularly
useful when the available data is restricted to a limited area of
state-space, which can lead parametric models to underestimate the
uncertainty in dynamics at population levels (states) which have not
been observed.  One reason for the relative absence of nonparametric
approaches in the natural resource management context may be the lack
of existing approaches for adapting the non-parametric Bayesian models
previously proposed [@Munch2005] to a decision-theoretic framework.
Adapting a non-parametric approach requires modification of existing
methods for decision theory.  We have illustrated how this might be
done for a classic stochastic dynamic programming problem, opening the
door for substantial further research into how these applications might
be improved.



<!-- On mechanistic models -->
<!-- { Seems kind of irrelevant }

The complexity of ecological interactions and a lack of data contribute
greatly to both of the problems. This concern is particularly acute
in case of ecological tipping points [@Barnosky2012; @Scheffer2001],
which arise from feedbacks common in ecological systems and can lead to a
sudden catastrophic transition to an undesirable state.  In general we do
not know where such tipping points are unless we have already observed
the transition, in which case it is frequently too late to respond.
Though there may exist early warning signals for certain kinds of these
transitions that are driven by slow changes [@Scheffer2009], we do not
know when, where, or how to apply them to the decision making context
more generally [@Boettiger2013]. 

--> 

----------------------------------------






Code to replicate the analysis, along with complete log of this research can be found at: [https://github.com/cboettig/nonparametric-bayes](https://github.com/cboettig/nonparametric-bayes/)


## Markov Chain Monte Carlo Analysis

```{r appendixplots, dependson="assemble-models"}
```


