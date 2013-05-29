% Approximate Dynamic Programming Introduction
% Carl Boettiger
% May 28 2013


# Why ADP

* ADP for the _curse(s) of dimensionality_

* An _optimizing simulator_: Better to have accurate detailed dynamics and an approximate solution, then optimal solution to approximate dynamics...

* A mathematical compression 

# Three curses of dimensionality

1. State space
1. Action space
1. Outcome space


# Step 0

Initialization

#### 0a. 
  Initialize $\bar V_t^0(S_t) $ for all states $S_t$

#### 0b. 
  Choose an initial state $S_0^1$

#### 0c. 
  Set n = 1




# Step 1

choose a random sample path $\omega^n$



# Step 2

For $t = 0, 1, 2 ... , T$ do:

#### 2a. 
  Solve: 
 
$$V_t(S_t) = \max_{x_t \in \chi_t} \left(C(S_t, x_t) + \gamma \sum_{s^{\prime} \in \mathcal{S}} \mathbb{P}(s^{\prime} | S_t^n, x_t) V_{t+1}^{n-1} (s^{\prime}) \right)$$

  and let $x^n_t$ be the value of $x_t$ that solves the maximization (argmax).  

####  2b. 
  
  Update $\bar V_t^{n-1}(S_t)$ using 


$$V_t^n(S_t) = \begin{cases} 
\hat v_t^n & S_t = S_t^n \\
\bar V_t^{n-1}(S_t) & \textrm{otherwise} 
\end{cases}$$

#### 2c. 
  
  Compute  $S^n_{t+1} = S^M(S_t^n, x^n_t, W_{t+1}(\omega^n))$



# Step 3

Let $n = n+1$.  If $n < N$ go to step 1.  


# As R Code

```r
for(n in 1:N){
  omega_n <- rlnorm(Tmax, 0, sigma_g)
  S_current <- S_0 
```
 explores faster when this is random


# As R Code


```r
for(t in 1:Tmax){
  s <- which.min(abs(S-S_current)) 
``` 

Where `s` is the index of the state we're considering

# As R Code

Find the action maximizing the value in the next step

$$V_t(S_t) = \max_{x_t \in \chi_t} \left(C(S_t, x_t) + \gamma \sum_{s^{\prime} \in \mathcal{S}} \mathbb{P}(s^{\prime} | S_t^n, x_t) V_{t+1}^{n-1} (s^{\prime}) \right)$$

```r   
values <- sapply(1:length(chi), function(x)
  C(S[s], chi[x]) + gamma * 
  sdp_matrix[[x]][s,] %*% V[,t])
hat <-  c(x_nt = which.max(values), 
          v_nt = max(values))
```

# As R Code

Update value V 

```r
V[hat["x_nt"], t] <- hat["v_nt"] 
```    


and advance the state in time along random path  

$$ S^n_{t+1} = S^M(S_t^n, x^n_t, W_{t+1}(\omega^n))$$

```r
S_current <- omega_n[t] * 
             f(S_current, chi[hat["x_nt"]], p)

  }
}
```



# What could go wrong?

We no longer have the loop-over-all-states problem, but we face several new or remaining issues:

1. We still require the use of the one-step transition matrix, with the equally troublesome sum over all states $\sum_{s^{\prime}\in S} \mathbb{P}(s^{\prime} | S_t^n, x_t)$.  
2. We only update the values of states we visit.  We still need a way to estimate the value of states we have not visited.  
3. Worse, we might not visit states that seem bad relative to states we have visited.


# Stochastic value function sampling 

Dealing with problem 1:
 
$$V_t(S_t) = \max_{x_t \in \chi_t} \left(C(S_t, x_t) + \gamma \sum_{\hat \omega \in \hat \Omega^n_{t+1}} p_t+1(\hat \omega) \bar V_{t+1}^{n-1} (S_{t+1}) \right)$$

```r
    V[hat["x_nt"], t] <- (1 - alpha) * 
                         V[hat["x_nt"],t] + 
                          alpha * hat["v_nt"] 
```




# Decreasing the state space

* Aggregation
* Continuous Value function approximations
* Using the post-decision state variable



# Initialization problem 

* We do not explore if we are too pessimistic about value of visiting other states.  Start optimisitc: AI's A* alogrithm (synchronous)
* Asynchronous updating -- randomly sampling starting variables
* RTDP (Real Time Dynamic Programming -- not necessarily what it sounds like) external rule determines which states we visit


# Learning



