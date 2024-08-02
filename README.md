# SAS_GENGAMMA_QUANTILE
Shows closed-form computation of Generalized Gamma quantile using the stable parametrization of PROC LIFEREG dist = Gamma

# Background
The Generalized Gamma (Gengamma) CDF is not available in SAS. There is a relationship between the original Gengamma distirbution and the "standard" Gamma distribution. One can leverage on this to relationship to draw a Gengamma CDF and its inverse (quantile function) with applications to sampling using inverse transform sampling as shown [here](https://blogs.sas.com/content/iml/2021/03/15/generalized-gamma-distribution.html#comment-556361).

# Open issue 
Often a different Gengamma parametrization is used which offers increased stability during optimization. Therefore, a "stable" Gengamma distribution is used in parametric survival regression (PROC LIFEREG dist = Gamma). Here, the relationship to the "standard" Gamma distribution is not directly applicable due to the different parametrization. Hence, an closed-form quantile function for the stable Gengamma is not avaialble. The SAS [documentation](https://documentation.sas.com/doc/en/pgmsascdc/9.4_3.4/statug/statug_lifereg_details12.htm) states that predictions (e.g., quantile calculations) for the stable Gengamma are performed numerically and no closed-form solution is given. Here PROC LIFEREG will provide predictions within performing of a regression analysis.

# Task
Sometimes however we might want to use an existing regression model to predict unavailable survival times, or we might want to use the Gengamma quantile function flexibly for simulations purposes. In all such cases we need to write down the quantile functions in closed-form or have a viable method to perform the quantile operation. A quantile function for the stable Gengamma distribution is therefore needed.

# Solution
One can show that to solve this problem the missing step is to find the relationship between the original and stable Gengamma distribution. While not documented everywhere, one can show with some simple algebra that the original and stable Gengamma parametrizations are related as following (note the terminological conflicts):

|Parameter    | Stable as in PROC LIFEREG     | Original (as in [Lawless](https://www.jstor.org/stable/1268326))  |
|-------------| ------------- | ------------- |
|  scale ($\mu$ in SAS)          | $\mu = X\beta$  | $a = exp(\mu - log(k)/b)$  |
|  free shape ($\delta$ in SAS)      |   $\delta$        | $k = \delta^{-2}$  |
|  shape ("scale" in SAS)          |  $\sigma$           | $b = \delta / \sigma$     | 

Where $\beta$ is the  vector of coefficient and $X$ the design matrix in PROC LIFEREG.

# Application
Given the above relationship one can show that the original Gengamma desity is equal to the "standard" Gamma, $(\frac{x}{a})^b \sim G(k, 1)$ with shape = $k$ and scale = 1. It follows that we can analytically calculate a stable Gengamma log quantile as 
```math
log(q_p) = log(a) + \frac{1}{b} \cdot log\left( G_{k,1}^{-1}(p) \right) 
```
replacing $a$, $k$, and $b$ with the above definitions where $G_{k,1}^{-1}(p)$ is the inverse "standard" Gamma with parameters $k, 1$ and typically the percentile $p = 0.5$.
