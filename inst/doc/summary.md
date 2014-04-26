The fishery manager's problem is often posed on a infinite time horizon as follows;

$$\max_{q_t \geq 0} \mathbb{E} \lbrace \sum_0^{\infty} \alpha^t \Pi(h_t, x_t) \rbrace$$

s.t.

\begin{align*}
x_t = z^g G(s_{t-1}) \\
s_t = x_t - h_t\\
m_t = z_t^m x_t\\
h_t = \min(x_t z_t^i q_t)
\end{align*}

Where $h$ is harvest, $x$ the stock size (we refer to the "fish stock", not a wallstreet stock ;-), which is subject to stochastic growth shocks $z^g_t$.  Growth of fish depends on the population you didn't harvest, $s_t$.  The stock size $x_t$ is only known with some uncertainty, hence the decision-maker only sees $m_t$ measurement shock $z^m_t$.  The decision-maker chooses a quota $q_t$ at each time step, which is also implemented with error as the realized harvest $h_t$, and $\alpha$ is the discount rate, and $\Pi$ the profit as a function of how many fish we harvest $h_t$.  So we can write Bellman recursion as 

$$J_t(m_t) = \max_{q_t \geq 0} \mathbb{E} \lbrace h_t + \alpha J_{t+1} (z_{t+1}^m z_{t+1}^g G(x_t = h_t))\rbrace$$

Typically this might be solved numerically by stochastic dynamic programming given the equation of motion / dynamical system $G$.  This requires making everything discrete of course, which I do just using a uniform grid but no doubt some more clever basis functions could be used.  

The problem becomes more challenging when $G$ must be estimated from the data first.  Typically this is done using a particular parametric equation and things are straightforward: e.g. we get posterior densities on the parameters of this equation and we integrate over them when calculating the expectation in the recursion equation.  

So my interest steps in where we don't even have a structural equation for $G$, we just have some time series of previous observed dynamics (i.e. measured stock-size $m$. along with perhaps a previous history of the action/policy $q$).  So instead we want to construct a model of the data from scratch, e.g. with a machine-learning approach.  We want to capture the uncertainty in our estimate of G, along with the intrinsic stochasticity/errors in the process, measurement, and implementation (the $z$ shocks above -- also of unknown parameteriziation.  

I have taken non-parametric Bayesian approach so far using Gaussian process regression, (ignoring the measurement and implementation error for the moment, which allows me to solve in the space of $x=m$ and $h=q$.)  This works reasonably well in the cases I have looked at so far; in particular the case of interest when the true dynamics for $x$ have a fold bifurcation (tipping point) in a previously unmeasured area of state space.  (e.g. some examples in my [notes](http://carlboettiger.info/2013/04/27/comparison-of-nonparametric-and-parametric-approaches.html)).  

I'd be curious to hear more about standard approaches to this case where we lack the underlying dynamical model.  I'd also be interested in knowing more about the dual-control problem (something we ecologists call active adaptive management) where the decision variable (harvest quota) is intentionally manipulated to learn about the system (e.g. analgous to the one-armed bandit problem).  I've posed the question in 1 dimension, but in reality the fish dynamics are multi-dimensional (e.g. different age classes of fish, different interacting species, spatial structure, etc).  Obviously that makes the curse of dimensionality that much worse.  
