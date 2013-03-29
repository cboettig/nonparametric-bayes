Non-parametric approaches to optimal policy are more robust
===========================================================

Carl Boettiger, Marc Mangel, Steve Munch

Abstract
========

Introduction
============

The problem of structural uncertainty in managing ecological systems.

-   Most management recommendations from the ecological literature /
    management policies based on (motivated by) parametric models.
    *Preference to pitch towards policy or theoretical literature?*

-   Background discussion on the importance/success of parametric
    modeling; (e.g. Levins (1966), etc., up through Geritz and Kisdi
    (2011))

-   Background on the concerns of structural uncertainty – we don’t have
    the right models. (Also: measurement uncertainty, parameter
    uncertainty, unobserved states, boundary conditions, etc.)

-   Hierarchical Bayesian approach has provided a natural way to address
    these from a statistical standpoint. Successes and challenges from
    the parametric route: e.g. Cressie et al. (2009).

-   Among pitfalls of these approaches: particularly difficult to apply
    in management context. Optimization-based (Decision-theoretic)
    approaches to ecological management need to be precisely
    parameterized about everything, even the uncertainty itself (Polasky
    et al. 2011) (discussion of Brozović and Schlenker (2011) SDP
    example in context of threshold system would perhaps be useful as
    well).

-   Management goals / decision-theoretic approaches need accurate
    prediction over relevant (short?) timescales more than accurate (but
    incomplete or noisy estimated) mechanisms.

-   Nonparametric (machine-learning?) approaches may offer the benefit
    of the hierarchical Bayesian approach without the practical and
    computational limitations of their parametric kin. Non-parametric
    models are flexible enough to take maximum advantage of the data
    available, while being appropriately ambiguous about the dynamics of
    a system in regions of parameter space that have been poorly or
    never sampled.

-   Nonparametric approaches are beginning to appear more frequently in
    ecological and conservation literature (Species distribution
    models/maxent), including the Gaussian process based approach used
    here (Munch et al. 2005). However, such approaches have yet to be
    applied to the decision-theoretic framework that could guide
    management decisions, where we expect them to excel for several
    reasons (1) the ability to make accurate forecasts by more closely
    approximating the underlying process where data is available (2)
    remaining appropriately ambiguous where data is not available (3)
    remaining computationally simple enough to avoid some pitfalls
    common to hierarchical parametric approaches.

Approach and Methods
====================

### (Brief) Discussion of GP inference

-   The use of Gaussian processes for inference in dynamical systems
    [introduced by @Kocijan2005]

-   Gaussian processes in ecological literature (Munch et al. 2005)

-   Our methodology (e.g. following Munch et al. (2005)).

### Discussion of the dynamic programming solution

*(More thorough, but general-audience targeted. Technical details and
code provided in appendices).*

The fishery management problem over an infinite time horizon can be
stated as:

\begin{align}
& \max_{ \{h_t\} \geq 0 } \mathbf{E} \lbrace \sum_0^\infty \delta^t \Pi(h_t) \rbrace \\
& \mathrm{s.t.}  \\
 & X_t = Z_t f\left(S_{t-1}\right) \\
 & S_t = X_t - h_t \\
 & X_t  \geq 0 
\end{align}
Where $\mathbf{E}$ is the expectation operator, $\delta$ the discount
rate, $\Pi(h_t)$ the profit expected from a harvest of $h_t$, and other
terms as in Eq. (1). For simplicity, we have assumed that profits depend
only on the chosen harvest; simplifying further we will usually consider
profits to be proportional to harvest, $\Pi(h_t) = h_t$.

Once the posterior Gaussian process (GP) has been estimated (e.g. see
Munch et al. 2005), it is necessary to adapt it in place of the
parametric equation for the stochastic dynamic programming (SDP)
solution (see Mangel and Clark 1988 for a detailed description of
parametric SDP methods) to the optimal policy. The essense of the idea
is straight forward – we will use the estimated GP in place of the
parametric growth function to determine the stochastic transition matrix
on which the SDP calculations are based.

The posterior Gaussian process is completely defined by the expected
value and covariance matrix at a defined set of training points. For
simplicty we will consider a these points to fall on a discrete, uniform
grid $x$ of 101 points from 0 to 15 (1.5 times the positive equilibrium
$K$). Again to keep things simple we will use this same grid
discritization for the parametric approach. Other options for choosing
the grid points, including collocation methods and functional basis
expansion (or even using Guassian processes in place of the discrete
optimization; an entirely different context in which GP can be used in
SDP, see [@Deisenroth2009]) could also be considered.

The transition matrix $\mathbf{F}$ is thus an 101 by 101 matrix for
which the ${i,j}$ entry gives the probability of transitioning into
state $x_i$ given that the system is in state $x_j$ in the previous
timestep. To generate the transition matrix based on the posterior GP,
we need only the expected values at each grid point and the
corresponding variances (the diagonal of the covariance matrix), as
shown in Figure 1. Given the mean at each gridpoint as the length 101
vector $E$ and variance $V$, the probability of transitioning from state
$x_i$ to state $x_j$ is simply
$\mathcal{N}\left(x_j | \mu = E_i, \sigma = \sqrt{V_i}\right)$, where
$\mathcal{N}$ is the Normal density at $x_j$ with mean $\mu$ and
variance $\sigma^2$. Strictly speaking, the transition probability
should be calculated by integrating the normal density over the bin of
width $\Delta$ centered at $x_j$. For a sufficiently fine grid that
$f(x_j) \approx f(x_j + \Delta)$, it is sufficient to calculate the
density at $x_j$ and then row-normalize the transition matrix.

Pseudocode for the determining the transtion matrix from the GP
---------------------------------------------------------------

~~~~ {.r}
for(h in h_grid)
  F_h = for(x_j in grid)
          for(i in 1:N) 
            dnorm(x_j, mu[i]-h, V[i])
~~~~

A transition matrix for each of the parametric models $f$ is calculated
using the log-normal density with mean $f(x_i)$ and log-variance as
estimated by maximum likelihood. From the discrete transition matrix we
may write down the Bellman recursion defining the the stochastic dynamic
programming iteration:

\begin{equation}
V_t(x_t) = \max_h \mathbf{E} \left( h_t + \delta V_{t+1}( Z_{t+1} f(x_t - h_t)) \right)
\end{equation}
where $V(x_t)$ is the value of being at state $x$ at time $t$, $h$ is
control (harvest level) chosen. Numerically, the maximization is
accomplished as follows. Consider the set of possible control values to
be the discrete 101 values corresponding the the grid of stock sizes.
Then for each $h_t$ there is a corresponding transition matrix
$\mathbf{F}_h$ determined as described above but with mean
$\mu = x_j - h_t$. Let $\vec{V_t}$ be the vector whose $i$th element
corresponds to the value of having stock $x_i$ at time $t$. Then let
$\Pi_h$ be the vector whose $i$th element indicates the profit from
harvesting at intensity $h_t$ given a population $x_i$ (e.g.
$\max(x_i, h_t)$ since one cannot harvest more fish then the current
population size). Then the Bellman recursion can be given in matrix form
as

$$V_{t} = \max_h \left( \Pi_{h_{t}} + \delta \mathbf{F}_h V_{t+1} \right)$$

where the sum is element by element and the expectation is computed by
the matrix multiplication $\mathbf{F} V_{t+1}$.

### Pseudocode for the Bellman iteration

~~~~ {.r}
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
~~~~

*Currently this shows the literal R code, should be adapted*

### Discussion on how we compare performance of policies

-   Replicate stochastic simulations
-   Sensitivity analysis (Figure 4).

Example in a bistable system
----------------------------

Concerns over the potential for tipping points in ecological dynamics
(Scheffer et al. 2001) highlight the dangers of uncertainty in
ecological management and pose a substantial challenge to existing
decision-theoretic approaches (Brozović and Schlenker 2011). To compare
the performance of nonparametric and parametric approaches in an example
that is easy to conceptualize, we will focus on a simple parametric
model for a single species (derived from fist principles by Allen et al.
2005) as our underlying “reality”.

\begin{align}
X_{t+1} &= Z_t f(S_t) \\
S_t &= X_t - h_t \\
f(S_t) &= e^{r \left(1 - \frac{S_t}{K}\right)\left(S_t - C\right)}
\end{align}
Where $Z_t$ is multiplicative noise function with mean 1, representing
stochastic growth. We will consider log-normal noise with shape
parameter $\sigma_g$. We start with an example in which the parameters
are $r =$ 2, $K =$ 10, $C =$ 5, and $\sigma_g =$ 0.05.

As a low-dimensional system completely described by three parameters,
this scenario should if anything be favorable to a parametric-based
approach. This model contains an Allee effect, or tipping point, below
which the population is not self-sustaining and shrinks to zero
(Courchamp, Berec, and Gascoigne 2008).

Sample training data
--------------------

Both parametric and nonparametric approaches will require training data
on which to base their model of the process. We generate the training
data under the model described in Eq 1 for 35 time steps, under a known
but not necessarily optimal sequence of harvest intensities, $h_t$. For
simplicity we imagine a fishery that started from zero harvest pressure
and has been gradually increasing the harvest.

<!--Should we include any emprical examples? -->
Using data simulated from a specified model rather than empirical data
permits the comparison against the true underlying dynamics, setting a
bar for the optimal performance possible.

(Motivation, alternatives, stationarity, examples without a stable node
(limit-cycle models), examples based on observations near a stable node
alone, and why that isn’t impossible).

Results
=======

### Discussion of maximum likelihood estimated models

We estimate two parametric models from the data using a maximum
likelihood approach. The first model is structurally identical to the
true model (Eq 1), differing only in that it’s parameters are estimated
from the observed data rather than given. The alternative model is the
Ricker model, which is structurally similar and commonly used in for
such data.

(MLE models will assume the noise is log-normal, which it is in the
simulation).

Which estimates a Ricker model with $r =$ 1.8501, $K =$ 9.8091, and the
Allen Allele model with $r =$ 2.8079, $K =$ 11.8235 and $C =$ 7.2159.

Figure 1:
---------

*Shows the inferred Gaussian Process compared to the true and parametric
models. Refer to the appendix for details on the GP posteriors, etc.*

![Graph of the inferred Gaussian process compared to the true process
and maximum-likelihood estimated process. Graph shows the expected value
for the function $f$ under each model. Two standard deviations from the
estimated Gaussian process covariance with (light grey) and without
(darker grey) measurement error are also shown. The training data is
also shown as black points. (The GP is conditioned on 0,0, shown as a
pseudo-data
point).](http://carlboettiger.info/assets/figures/2013-02-12-20-21-16-a15dc1e581-figure/gp_plot.pdf)

Figure 2:
---------

*The take-home message, showing that the GP is closest to the optimal
strategy, while the parametric methods are less accurate. Visualizing
the policy may be more useful for the technical reader, the general
audience may prefer Figure 3 showing all replicates of the population
collapse under the parametric model and not under the GP.*

![The steady-state optimal policy (infinite boundary) calculated under
each model. Policies are shown in terms of target escapement, $S_t$, as
under models such as this a constant escapement policy is expected to be
optimal (Reed
1979).](http://carlboettiger.info/assets/figures/2013-02-12-20-21-26-a15dc1e581-figure/policies_plot.pdf)

Figure 3:
---------

*Figure 3 is a less abstract and more visceral visualization of the
take-home message, with the structurally inaccurate model leading
universally to a collapse of the fishery and very few profits, while the
Gaussian process performs nearly optimally. The parametric approach even
with the correct underlying structure does not perform optimally,
choosing in this case to under-fish (may need to show harvest dynamics
since that is not clear from the figure! Also isn’t general, sometimes
does optimally, sometimes over-fishes. Perhaps need to show more
examples.) May need to show profits too?*

![Gaussian process inference outperforms parametric estimates. Shown are
100 replicate simulations of the stock dynamics (eq 1) under the
policies derived from each of the estimated models, as well as the
policy based on the exact underlying
model.](http://carlboettiger.info/assets/figures/2013-02-12-20-21-33-a15dc1e581-figure/sim_plot.pdf)

Figure 4:
---------

*Shows the sensitivity analysis. A histogram of distribution of yield
over stochastic realizations, showing that the qualitative results do
not depend on the stochastic realization of the training data here, or
on the parameters of the underlying model, though quantitative
differences are visible.*

Discussion / Conclusion
=======================

-   Non-parametric methods have received far too little attention in
    ecological modeling efforts that are aimed at improved conservation
    planning and decision making support.

-   Importance of non-parametric approaches in conservation planning /
    resource management / decision theory.

-   Decision-theoretic tools such as optimal control calculations rely
    on robust *forecasting* more strongly than they rely on accurate
    *mechanistic* relationships.

-   Adapting a non-parametric approach requires modification of existing
    methods for decision theory. We have illustrated how this might be
    done in the context of stochastic dynamic programming, opening the
    door for substantial further research into how these applications
    might be improved.

-   Anticipate improved relative performance in higher dimensional
    examples

-   Discuss constant escapement in model, in policies.

-   Limitations of this comparison: Are the maximum-likelihood solutions
    a straw man?

-   Discussion of alternative related approaches: POMDP/MOMDP,

Future directions
-----------------

-   Multiple species \* Online learning \* Multiple step-ahead
    predictions \* Explicitly accomidating additional uncertainties \*
    Improving inference of optimal policy from the GP

`{r echo=FALSE, results="asis"} #bibliography("html")`

Appendix / Supplementary Materials
==================================

MCMC posterior distributions and convergence analysis
-----------------------------------------------------

![Histogram of posterior distributions for the estimated Gaussian
Process shown in Figure 1. Prior distributions
overlaid.](http://carlboettiger.info/assets/figures/2013-02-12-20-21-38-a15dc1e581-figure/posteriors.pdf)

@Gramacy2005

Tables of nuisance parameters, sensitivity analysis
---------------------------------------------------

### List of hyper-parameters, prior distributions and their parameters

Reproducible code, “Research Compendium”
----------------------------------------

Allen, Linda J. S., Jesse F. Fagan, Göran Högnäs, and Henrik Fagerholm.
2005. “Population extinction in discrete-time stochastic population
models with an Allee effect.” *Journal of Difference Equations and
Applications* 11 (4-5) (apr): 273–293. doi:10.1080/10236190412331335373.

Brozović, Nicholas, and Wolfram Schlenker. 2011. “Optimal management of
an ecosystem with an unknown threshold.” *Ecological Economics* (jan):
1–14. doi:10.1016/j.ecolecon.2010.10.001.

Courchamp, Franck, Ludek Berec, and Joanna Gascoigne. 2008. *Allee
Effects in Ecology and Conservation*. Oxford University Press, USA.

Cressie, Noel, Catherine a Calder, James S. Clark, Jay M. Ver Hoef, and
Christopher K. Wikle. 2009. “Accounting for uncertainty in ecological
analysis: the strengths and limitations of hierarchical statistical
modeling.” *Ecological Applications* 19 (3) (apr): 553–70.

Geritz, Stefan a. H., and Eva Kisdi. 2011. “Mathematical ecology: why
mechanistic models?” *Journal of mathematical biology* (dec).
doi:10.1007/s00285-011-0496-3.

Levins, Richard. 1966. “The strategy of model building in population
biology.” *American Scientist* 54 (4): 421–431.

Mangel, Marc, and Colin W. Clark. 1988. *Dynamic Modeling in Behavioral
Ecology*. Princeton: Princeton University Press.

Munch, Stephan B., Melissa L. Snover, George M. Watters, and Marc
Mangel. 2005. “A unified treatment of top-down and bottom-up control of
reproduction in populations.” *Ecology Letters* 8 (7) (may): 691–695.
doi:10.1111/j.1461-0248.2005.00766.x.

Polasky, Stephen, Stephen R. Carpenter, Carl Folke, and Bonnie Keeler.
2011. “Decision-making under great uncertainty: environmental management
in an era of global change.” *Trends in ecology & evolution* (may): 1–7.
doi:10.1016/j.tree.2011.04.007.

Reed, William J. 1979. “Optimal escapement levels in stochastic and
deterministic harvesting models.” *Journal of Environmental Economics
and Management* 6 (4) (dec): 350–363. doi:10.1016/0095-0696(79)90014-7.

Scheffer, Marten, Stephen R. Carpenter, J. A. Foley, C. Folke, and B.
Walker. 2001. “Catastrophic shifts in ecosystems.” *Nature* 413 (6856)
(oct): 591–6. doi:10.1038/35098000.
