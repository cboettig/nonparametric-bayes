---
title: "Avoiding tipping points in the management of ecological systems: a non-parametric Bayesian approach"
author: 
  - name: "Carl Boettiger"
    affiliation: cstar
    email: cboettig@gmail.com
    footnote: Corresponding author
  - name: Marc Mangel
    affiliation: cstar
  - name: Stephan Munch
    affiliation: noaa
address: 
  - code: cstar
    address: | 
      Center for Stock Assessment Research, 
      Department of Applied Math and Statistics, 
      University of California, Mail Stop SOE-2,
      Santa Cruz, CA 95064, USA
  - code: noaa
    address: | 
      Southwest Fisheries Science Center, 
      National Oceanic and Atmospheric Administration, 
      110 Shaffer Road, Santa Cruz, CA 95060, USA

abstract: |
          Model uncertainty and limited data coverage are fundamental challenges to
          robust ecosystem management.  These challenges are acutely highlighted
          by concerns that many ecological systems may contain tipping points.
          Before a collapse, we do not know where the tipping points lie, if the
          exist at all.  Hence, we know neither a complete model of the system
          dynamics nor do we have access to data in some large region of state-space
          where such a tipping point might exist.  These two sources of uncertainty
          frustrate state-of-the-art parametric approaches to decision theory
          and optimal control.  I will illustrate how a non-parametric approach
          using a Gaussian Process prior provides a more flexible representation
          of this inherent uncertainty.  Consequently, we can adapt the Gaussian
          Process prior to a stochastic dynamic programming framework in order to
          make robust management predictions under both model and uncertainty and
          limited data.
keywords:
  - Bayesian
  - Structural Uncertainty
  - Nonparametric 
  - Optimal Control
  - Decision Theory
  - Gaussian Processes
  - Fisheries Management
---

<!-- Run R code for analysis, to be called by figures -->



















































































































Introduction
============

Decision making under uncertainty is a ubiquitous challenge of
natural resource management and conservation.  Ecological dynamics are
frequently complex and difficult to measure, making uncertainty in our
understanding a prediction a persistent challenge to effective management.
Decision-theoretic approaches provide a framework to determine the
best sequence of actions in face of uncertainty, but only when that
uncertainty can be meaningfully quantified [@Fischer2009].  The sudden
collapse of fisheries and other ecosystems has increasingly emphasized
the difficulties of formulating even qualitatively correct models of
the underlying processes.


We develop a novel approach to address these concerns in the context of fisheries; though
the underlying challenges and methods are germane to many other
conservation and resource management problems.  The economic value and
ecological concern have made marine fisheries the crucible for much
of the founding work [@Gordon1954; @Reed1979; @May1979; @Ludwig1982]
in managing ecosystems under uncertainty.  Global trends [@Worm2006]
and controversy [@Hilborn2007; @Worm2009] have made understanding these
challenges all the more pressing.

<!-- Uncertainty outside the data without the correct model has not been handled. -->
<!-- __We don't have the model__ -->

Uncertainty enters the decision-making process at many levels:
intrinsic stochasticity in biological processes, measurements, and
implementation of policy [_e.g._ @Reed1979; @Clark1986; @Roughgarden1996;
@Sethi2005], parametric uncertainty [_e.g._ @Ludwig1982; @Hilborn1997;
@McAllister1998; @Schapaugh2013], and model or structural uncertainty
[_e.g._ @Williams2001; @Cressie2009;  @Athanassoglou2012].  Of these,
structural uncertainty incorporates the least a priori knowledge or
assumptions and is generally the hardest to quantify. Typical approaches
assume a weak notion of model uncertainty in which the correct model
(or reasonable approximation) of the dynamics must be identified from
among a handful of alternative models.  Here we consider an approach
that addresses uncertainty at each of these levels without assuming the
dynamics follow a particular (i.e. parametric) structure.



_Cut the next three paragraphs, since they are covered more concisely in the above paragraph?_

### Process, measurement, and implementation error

Resource management and conservation planning seek to determine the
optimal set of feasible actions to maximize the value of some objectives
(e.g @Halpern2013).  Process error, measurement error, implementation
error [@Reed1979, @Clark1986, @Roughgarden1996, @Sethi2005].  These
sources of stochasticity in turn mean that model parameters can only
be estimated approximately, requiring parametric uncertainty also be
considered [@Ludwig1982].

### Parametric uncertainty

As the parameter values for these models must be estimated from limited
data, there will always be some uncertainty associated with these values.
This uncertainty further compounds the intrinsic variability introduced
by demographic or environmental noise.  The degree of uncertainty in
the parameter values can be inferred from the data and reflected in the
estimates of the transition probabilities [@Ludwig1982; @Mangel1988;
@Hilborn1997; @Schapaugh2013].

### Structural (model) uncertainty

Estimates of parameter uncertainty are only as good as the parametric
models themselves.  Often we do not understand the system dynamics well
enough to know if a model provides a good approximation over the relevant
range of states and timescales (criteria that we loosely refer to as
defining the "right" or "true" model.)  So called structural or model
uncertainty is a more difficult problem than parametric uncertainty.
Typical solutions involve either model choice, model averaging, or
introducing yet greater model complexity of which others may be special
cases (model averaging being one such way to construct such a model)
[@Williams2001; @Athanassoglou2012; @Cressie2009].  Even setting aside
other computational and statistical concerns (e.g. [@Cressie2009]), these
approaches do not address our second concern - representing uncertainty
outside the observed data range.


<!-- __We don't have the data where we need it__ -->
<!-- What do we call this?  Extrapolation uncertainty?  Pathological Uncertainty? -->


Model uncertainty is particularly insidious when model predictions must
be made outside of the range of data on which the model was estimated.
This extrapolation uncertainty is felt most keenly in decision-theoretic
(or optimal control) applications, since (a) exploring the potential
action space typically involves considering actions that may move the
system outside the range of observed behavior, and (b) decision-theoretic
alogrithms rely not only on reasonable estimates of the expected
outcomes, but depend on the weights given to all possible outcomes
[_e.g._ @Weitzman2013].  If we are observing the fluctuations of a
given fish stock over many years under a fixed harvesting pressure,
we might develop and test a model that could reasonably predict the
frequency of a deviation of a given size, even when such a deviation has
not been previously observed. Yet such predictions are far less reliable
when extrapolated to a harvest pressure that has not yet been observed.
Thus, model uncertainty can be particularly challenging in the management
and decision-making context.


This difficult position of having neither the true model nor data that
covers the full range of possible states is unfortunately the rule more
than the exception. The potential concern of tipping points in ecological
dynamics [@Scheffer2001; @Polasky2011] reflects these concerns -- as
either knowledge of the true model or more complete sampling of the
state space would make it easy to identify if a tipping point existed.
If we do not know but cannot rule out such a possibility, then we
face decision-making under this dual challege of model uncertainty and
incomplete data coverage.


These dual concerns pose a substantial challenge to existing decision-theoretic
approaches [@Brozovic2011].  Because intervention is often too late 
after a tipping point has been crossed (but see @Hughes2013), management
is most often concerned with avoiding potentially catastrophic tipping
points before any data is available at or following a transition that
would more clearly reveal these regime shift dynamics [e.g. @Bestelmeyer2012].

<!--
_Map_
-->

Here we illustrate how a stochastic dynamic programming (SDP) algorithm
[@Mangel1988; @Marescot2013] can be driven by the
predictions from a Bayesian non-parametric (BNP) approach [@Munch2005a].
This provides two distinct advantages compared with contemporary
approaches.  First, using a BNP sidesteps the need for an accurate
model-based description of the system dynamics.  Second, the BNP can
better reflect uncertainty that arises when extrapolating a model outside
of the data on which it was fit.  We illustrate that when the correct
model is not known, this latter feature is crucial to providing a robust
decision-theoretic approach in face of substantial structural uncertainty.





<!-- necessary? --> <!-- should reflect the dual problem of extrapolation and model uncertainty better -->

This paper represents the first time the SDP decision-making framework has
been used without an a priori model of the underlying dynamics through
the use of the BNP approach.  In contrast to parametric models which
can only reflect uncertainty in parameter estimates, the BNP approach
provides a more state-space dependent representation of uncertainty.
This permits a much greater uncertainty far from the observed data than
near the observed data.  These features allow the BNP-SDP approach to
find robust management solutions in face of limited data and without
knowledge of the correct model structure.



<!--
-  _Note on "not magic"_: honest uncertainty + SDP
--> 

The idea that any approach can perform well without either having to
know the model or have particularly good data should immediately draw
suspicion.  The reader must bear in mind that the strength of our approach
comes not from black-box predictive power from such limited information,
but rather, by providing a more honest expression of uncertainty outside
the observed data without sacrificing the predictive capacity near the
observed data. By coupling this more accurate description of what is known
and unknown to the decision-making under uncertainty framework provided
by stochastic dynamic programming, we are able to obtain more robust
management policies than with common parametric modeling approaches.



<!--move to  much later *Or do we need this at all?* 
- _Note on comparing models_ (via value function rather than by "fit"). 
-->

The nature of decision-making problems provides a convenient way to compare 
models.  Rather than compare models in terms of best fit to data or fret over
the appropriate penalty for model complexity, model performance is defined 
in the concrete terms of the decision-maker's objective function, which we
will take as given. (Much argument can be made over the 'correct' objective
function, e.g. how to account for the social value of fish left in the sea
vs. the commercial value of fish harvested; see @Halpern2013 for further 
discussion of this issue.  Alternatively, we can always compare model performance
across multiple potential objective functions.)  The decision-maker
does not necessarily need a model that provides the best mechanistic understanding
or the best long-term outcome, but rather the one that best estimates the 
probabilities of being in different states as a result of the possible actions. 



## Background on the Gaussian Process

<!-- 
- Background on non-parametric modeling.  
-->


Addressing the difficulty posed by extrapolation without knowing the
true model requires a nonparametric approach to model fitting: one
that does not assume a fixed structure but rather depends on the size
of the data (e.g. non-parametric regression or a Dirichlet process).
This established terminology is nevertheless unfortunate, as (a)
this approach still involves the estimation of parameters, and (b),
Statisticians use non-parametric to mean both this property (structure
is not fixed by the parameters) and an entirely different (and probably
more familiar) case in which the model does not assume any distribution
(e.g. non-parametric bootstrap, order statistics).  Some literature
thus uses the term semi-parametric, which merely adds ambiguity to
the confusion.

This non-parametric property -- having a structure explicitly dependent on
the data -- is precisely the property that makes this approach attractive
in face of the limited data sampling challenges discussed above.
Having fit a parametric model to some data, the model is completely
described by the values (or posterior distributions) of it's parameters.
The non-parametric model is not captured by its parameter values or
distributions alone. Either the model scales with the complexity of the data
on which it is estimated (e.g. nonparametric hierarchical approaches such
as the Dirchlet process) or the data points become themselves part of the model
specification, as in the nonparametric regression used here.

<!--
- Definition
- Previous application
--> 

The use of Gaussian process (GP) regression (or "kriging" in the geospatial
literature) to formulate a predictive model is relatively new in the
context of modeling dynamical systems [@Kocijan2005], and was first introduced
in the context ecological modeling and fisheries management in @Munch2005.
An accessible and thorough introduction to the formulation and use of
GPs can be found in @Rasmussen2006.

<!--
- Why it is particularly suited to these two problems
- (Why this is a novel application thereof)
--> 

<!-- Useless 
The essence of the GP approach can be captured in the following
thought experiment: An exhaustive parametric approach to the challenge
of structural uncertainty might proceed by writing down all possible
functional forms for the underlying dynamical system with all possible
parameter values for each form, and then consider searching over this
huge space to select the most likely model and parameters; or using a
Bayesian approach, assign priors to each of these possible models and
infer the posterior distribution of possible models. The GP approach
can be thought of as a computationally efficient approximation to this
approach. GPs represent a large class of models that can be though
of as capturing or reasonably approximating the set of models in this
collection.  By modeling at the level of the process, rather than the
level of parametric equation, we can more concisely capture the possible
behavior of these curves.  In place of a parametric model of the dynamical
system, the GP approach postulates a prior distribution of (n-dimensional)
curves that can be though of as approximations to a range of possible
(parametric) models that might describe the data. The GP allows us to
consider probabilities on a large set of possible curves simultaneously.
-->

The posterior distribution for the hyper-parameters of the Gaussian 
process model are estimated by Metropolis-Hastings algorithm, again with
details and code provided in the Appendix.  @Rasmussen2006 provides
an excellent general introduction to Gaussian Processes and @Munch2005a 
first discusses their application in the context of population dynamics
models such as fisheries stock-recruitment relationships.


Approach and Methods
====================



### Statement of the optimal control problem

To illustrate the application of the BNP-SDP approach and compare to the
predictions of the alternative parametric models we focus on the classical
problem of selecting the appropriate harvest level given an observation
of the stock size in the previous year  [@Reed1979; @Ludwig1982; @Mangel1988].
Given this observation and the model (together with the parameter uncertainty) of the stock
recruitment process, the manager seeks to maximize the value of the 
fishery over a fixed time interval of 50 years at a discount rate
of 0.01. The value function (profits) at time $t$ depends on the true stock size
$x_t$ and the chosen harvest level $h_t$.  For simplicity we assume profit 
is simply proportional in the realized harvest (only
enforcing the restriction that harvest can not exceed available stock).  


### Parametric models

<!-- 
- Statement of the models
--> 

We consider three candidate parametric models of the stock-recruitment
dynamics: The Ricker model, the Allen model [@Allen2005], the Myers
model [@Myers1995]. The familiar Ricker model involves two parameters, corresponding
to a growth rate and a carrying capacity, and cannot support alternative
stable state dynamics (though as growth rate increases it exhibits a
periodic attractor that proceeds through period-doubling into chaos. We
will generally focus on dynamics below the chaotic threshold for
the purposes of this analysis.) The Allen model resembles the Ricker
dynamics with an added Allee effect parameter [@Courchamp2008], below
which the population cannot persist.  The Myers model also has three
parameters and contains an Allee threshold, but has compensatory rather
than over-compensatory density dependence (resembling a Beverton-Holt
curve rather than a Ricker curve at high densities.)  
We assume multiplicative log-normal noise perturbs the growth predicted 
by the each of the deterministic model skeletons described above. This 
introduces one additional parameter $\sigma$ that must be estimated by each
model. 

<!-- equations just in appendix? -->

As we simulate training data from the Allen model, we will
refer to this as the structurally correct model.  The Ricker model is
thus a reasonable approximation of these dynamics far from the Allee
threshold (but lacks threshold dynamics), while the Myers model shares
the essential feature of a threshold but differs in the structure. Thus
we have three potential parametric models of the stock dynamics.

We introduce parametric uncertainty by first estimating each of the
candidate models from data on unexploited stock dynamics following
some perturbation (non-equilibrium initial condition) over several time
steps. This training data could be generated in several different ways
(such as known variable exploitation rates, etc.), as long as it reflects
the dynamics in some limited region of state space without impacting
the problem.  We consider a period of 40 years of training data: long
enough that the estimates are not dependent on the particular realization,
while longer times are not likely to provide substantial improvement
(i.e. the results are not sensitive to this interval).  Each of the models
(described below) is fit to the same training data, as shown in Figure 1.


<!-- Bayesian inference of parametric models --> 

We infer posterior distributions for the parameters of each model
in a Bayesian context using Gibbs sampling (implemented in R [@RTeam]
using jags, [@R2jags]).  We choose uninformative uniform priors for all
parameters (See Appendix, Figures and tables, and the R code
provided). One-step-ahead predictions of these model fits are shown in
Figure 1.  While alternative approaches to the estimation of the posteriors
(such as integrating out the rate parameter $r$ analytically and then
performing a grid search over the remaining parameter space), the 
approach of using a standard Gibbs sampler routine is both more general
and representitive of common practice in estimating posteriors for such models. 
Each sampling is tested for Gelman-Rubin convergence and results are robust
to longer runs.  

<!-- SDP via parametric models -->

An optimal policy function is then inferred through stochastic dynamic
programming for each model given the posterior distributions of the
parameter estimates.  This policy maximizes the expectation of the value
function integrated over the parameter uncertainty. (code implementing
this algorithm provided in the Appendix).


### The Gaussian Process model

<!-- Statement of model -->  

<!-- FIXME ... say much more about the GP, as per Munch 2005a   -->
<!-- See much longer discusion above on "Background on GP".  Move some of that to here?  -->


We also estimate a simple Gaussian Process defined by a radial basis
function kernel of two parameters: $\ell$, which gives the characteristic
length-scale over which correlation between two points in state-space
decays, and $\sigma$, which gives the scale of the process noise by which
observations $Y_{t+1}$ may differ from their predicted values $X_{t+1}$
given an observation of the previous state, $X_t$. @Munch2005a gives an
accessible introduction to the use of Gaussian Processes in providing a
Bayesian nonparametric description of the stock-recruitment relationship.

<!-- Inference of the BNP model --> 


We use a Metropolis-Hastings Markov Chain Monte Carlo to infer posterior
distributions of the two parameters of the GP (Figure S13, code in
appendix), under weakly informative Gaussian priors (see parameters in
table S5). As the posterior distributions differ substantially from the
priors (Figure S13), we can be assured that most of the information in
the posterior comes from the data rather than the prior belief.

<!--
In contrast to the parametric models, this posterior distribution is still
conditional on the training data. As such, the uncertainty near the
observed data.
-->

<!-- SDP via the model --> 

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

<!-- 
### Figure 1: Fitted Models

- All models fit the data quite well
- Information criteria would pick the simple, incorrect model.
--> 








































