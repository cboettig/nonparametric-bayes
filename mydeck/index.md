---
title       : Approximate dynamic programming
subtitle    : 
author      : Carl Boettiger
job         : UC Santa Cruz
framework   : io2012        # {io2012, html5slides, shower, dzslides, ...}
highlighter : highlight.js  # {highlight.js, prettify, highlight}
hitheme     : tomorrow      # 
widgets     : mathjax       # {mathjax, quiz, bootstrap}
mode        : selfcontained # {standalone, draft}
---



## Three curses of dimensionality

1. State space
1. Action space
1. Outcome space

---


Solve: 
 
$$V_t(S_t) = \max_{x_t \in \chi_t} \left(C(S_t, x_t) + 
            \gamma \sum_{s^{\prime} \in \mathcal{S}} \mathbb{P}(s^{\prime} | S_t^n, x_t) V_{t+1}^{n-1} s^{\prime} \right)$$

  and let $x^n_t$ be the value of $x_t$ that solves the maximization (argmax).  

