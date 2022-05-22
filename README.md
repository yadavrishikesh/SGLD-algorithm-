### README files

<a rel="license" href="http://creativecommons.org/licenses/by/4.0/"><img alt="Creative Commons License" style="border-width:0" src="https://i.creativecommons.org/l/by/4.0/88x31.png" /></a><br />This work is licensed under a <a rel="license" href="http://creativecommons.org/licenses/by/4.0/">Creative Commons Attribution 4.0 International License</a>.

# Reference Paper: R., Huser, R. and Opitz, T. (2022) A flexible Bayesian hierarchical modeling framework for spatially dependent peaks-over-threshold data. Spatial Statistics 

This set of R codes simulate and do the inference of the the product mixtture models detailed below: 
# PARAMETERS for Gamma-Gamma model: Y(s)= alpha(s) X1(s) * X2 * X3(s), where: gamma(x) is the gamma function
# X1(s): i.i.d Weibull(1/beta1,1/gamma(1+beta1)), with shape 1/beta1 and scale 1/gamma(1+beta1)
# X2 spatially constant values (varying with respect to time only) with marginal Weibull(1/beta2,1/gamma(1+beta2))
# X3(s): marginally Inv-Gamma(beta3,beta3-1), with rate beta3-1 >0 and shape beta3 > 1, dependence structure governed by Gaussian copula with exponential correlation function exp(-h/rho), h>=0, with range rho>0
# alpha(s): scale parameter of the model with exp(gamma0 + gamma1 * Z1(s) + gamma2 * Z2(s) + gamma3 * Z3(s)), where Zi(s), i=1,2,3 are the spatial covariates 



