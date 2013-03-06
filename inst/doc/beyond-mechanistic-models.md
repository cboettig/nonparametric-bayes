Mechanistic models have long been the gold standard of theoretical modeling in ecology (e.g. see @Geritz2012 or @Cuddington2013).  Only by understanding the processes involved can we make reliable long term predictions and build an knowledge of cause and effect that guides the hypotheses we make, the data we collect, and the management decisions we make.  Process-based models, whether expressed in the language of mathematics or English, identify the connection between mosquitoes and the spread of malaria, or greenhouses gases and climate change, guiding our approach to understand and manage these threats. Despite this central importance, I argue here that ecologists would do well to give greater attention to the role non-mechanistic models can play in ecological management and decision making.  The value of these approaches is greatest in a context where decisions are made over short time horizons, and updated as new data becomes available. 



## The alternative to process-based models: pattern-based modeling

Before we proceed further it would be useful to define our terms and lay some greater context for the discussion.  In this essay, as we examine the advantages and disadvantages of mechanistic modeling in light of the alternatives, so it will be useful to better define what the alternatives are.  The alternative to mechanistic model or process-based description of the world has always been a correlative, or pattern-based one. Historically pattern-based modelling meant regression (usually linear regression) and resided in methodology rooted in statistics departments and was the primary focus of empericial ecologists, while mechanistic modeling meant dynamical systems (ODEs, PDEs, SDEs and their discrete kin), residing in methodology from the mathematics department and was the primary focus of theoretical ecologists.  

I believe these historical divides continue to color much of the literature today, with theorists more skeptical than empericists of statistical methodology and vice versa.  Meanwhile, the ground underneath has shifted.  Computational power and heirarchical statistical methods (particularly approximate likelihood or simulation techniques such as particle filters and approximate Bayesian computing) have brought the inference richer dynamical systems into their fold, while pattern-based modeling has spawned an entirely new approach under the banner of machine learning that has devided the statistics community (See @Bremman2001, "The Two Cultures").  Before these changes, pattern-matching approaches could be critiqued on the grounds of their simplicity (not everything is linear) while mechanisitic approaches could be critiqued on their tenous connections to data -- model parameter values such as death rates, birth rates, etc would be estimated in advance and then stuck into the model, rather than estimated in the context of the model itself.  Today both of these critiques 



A simple example will illustrate many of these points.  

Mechanistic modeling emphasizes the importance of capturing the correct gross properties of a system over the tracking minute fluctuations.  For instance, in selecting and parameterizing model of a population of conservation concern, we may be most interested in getting the long-term behavior correct -- such as identifying if the dynamics support persistence of the population -- rather than worrying how well they reflect the year-to-year fluctuations.  We would certainly have good reason to prefer such a model over alternatives which are irreconcilable to the most basic biological processes, such as unbounded growth, or growth curves that do not pass through the origin in a closed system. So it may come as a surprise to realize that such obviously wrong models can perform as well or even better than reasonable mechanistic models in guiding ecological management and decision making.  

An example of such a comparison is shown in Figure 1, in which the optimal management decision is determined by stochastic dynamic programming algorithm using the biologically (a) plausible and (b) the implausible model of the population growth.  

Three elements contribute to the biologically implausible phenomenological model performing better in this context:   

### 1. Relevant state space

The dynamics occur over a range of state-space in which the biologically implausible model performs as well or better than the more plausible model.  This is the most obvious and immediate reason why the shortcomings that appear to make the model implausible do not make it useless.  This effect is enhanced in our example because the management actions help drive the system towards rather than away from this region of state-space.  

### 2. Predictive accuracy

This application relies only on the predictive accuracy of the model, not an interpretation of the parameters.  Predictive accuracy is not the goal of all modeling, as ecologists have been observing for as long as they made models (perhaps none more memorably than @Levins1969).  

Mechanistic modeling is at its most powerful not when it is used to directly forecast future states but when it provides an understanding of how to approach a problem.  SIR-type models from the epidemiological literature are a good example.  While the simplest SIR models have little predictive power over the outbreak intensity or timing at a particular location, they provide a powerful understanding of the spread of an infection in terms of a single, biologically meaningful parameter: $R_0$, the basic reproductive number.  From the model, it becomes clear that management need not vaccinate every member of the population to stop the spread, but rather it suffices to vaccinate a sufficient fraction of the population to reduce $R_0$ below 1. 

### 3. Time scale for new data 

In the sequential decision making problem we considered, we are presented with new data after each action.  The relevant timescale for the prediction is thus not the long-term dynamics, which would be wildly divergent, but the dynamics over this much shorter interval.  While ecologists may be hesitant to base continual management on a model with obviously inaccurate long-term behavior, engineers tend to consider the problem in frequency space and gravitate to the opposite position -- a good control model should prioritize high-frequency accuracy over low frequency accuracy.  The differences in intuition may arise from the timescales at which each profession can typically adjust their control variables -- much faster for a control system of a chemical plant than state policy for a natural resource.  Still, the lesson is clear: when facing repeated management decisions over a short timescale, such as setting annual harvests of a fishery, it may be more valuable to use a machine learning algorithm that makes accurate year-ahead predictions that capture some of the high-frequency fluctuations that appear only as noise in a mechanistic model of the long-term population dynamics.  






## Advantages of non-mechanistic models

1. Can better use all of the data available.  

Incoprating various sources of information into mechanistic models can be an immensely difficult due to the increased complexity involved.  Only thanks to tandem advances in increasing computational power and hierarchical statistical methodology have we been able to tackle such intuitively important complexity (and the potentially new available data that accompanies it) such as spatial distribution, heterogeneties of space, time, and individuals, to shift to ecosystem-based approaches from single-species based approaches.  Machine learning algorithms 

2. Better express ambiguity in regions where data is not available

One of the greatest stengths of mechanistic models is their greatest weakness as well.  


## Better the devil you know?  The dangers of phenomenological approaches

While non-mechanistic models from the machine learning literature have the 

Few examples illustrate the danger of phenomenological modeling more broadly than the sub-prime mortgage crisis of 2008.  Simplifying greatly, the unquenchable demand for mortgage backed securities had shifted the underlying dynamics on which the phenomenological models were based faster than the models could learn about the shift.  Historical data indicated that mortgage-backed securities involved little risk.  


Mechanistic models are not free from this danger either. In fact they can be more suseptible to them for the very same reasons.  A mechanistic model has the ability to more convincingly produce out-of-sample predictions or extrapolation.  This may allow us to determine that a virus will spread or an endagered species will go extinct even before we have observed 


Biased estimation.  


Mechanistic models have their greatest strengths in a data-limited world -- where knowledge of plausible mechanisms can greatly decrease the risk 




Parametric models are frequently used phenomenologically, rather than being derived by plausible underlying mechanisms.  @Geritz2012 condemn this practice on the basis of the common misinterpretations that result.  These mistakes arise from ignoring or misinterpreting the biological meaning of the model parameters. The use of delay equations is a frequent offender -- for instance, one is hard-pressed to describe a mechanistic model consistent with Hutchinson's 1948 delay-logistic.  Assuming either scramble competition or interference competition gives rise to two alternative delay equations, both structurally different from the original.  

Three elements contribute to this.  

phenomenological accuracy over 


$$ $$


A process-based understanding helps us formulate hypotheses, guide experimental design and data collection, predict outcomes, and influence management decisions.  The alternative to a processed-based model is  a correlative one. Statical models machine learning.  




> since the late 1980s there has been some consensus amongst ecologists that management decisions are best guided by models which are grounded in ecological theory, and which strike a balance between too much or too little detail describing the relevant processes (e.g., DeAngelis 1988, Starfield 1997, Jackson et al. 2000, Carpenter 2003, Nelson et al. 2008) 

machine learning algorithms 

Prediction without 

Example of Statistical model vs process-based model.  




@Geritz2012 @Cuddington2013 


Parametric doesn't mean process-based/mechanistic



Managing by rule-of-thumb

Managing by algorithm


