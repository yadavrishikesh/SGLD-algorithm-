### README files

<a rel="license" href="http://creativecommons.org/licenses/by/4.0/"><img alt="Creative Commons License" style="border-width:0" src="https://i.creativecommons.org/l/by/4.0/88x31.png" /></a><br />This work is licensed under a <a rel="license" href="http://creativecommons.org/licenses/by/4.0/">Creative Commons Attribution 4.0 International License</a>.

### For more details about the R code, contact rishikesh.yadav@kaust.edu.sa and see the below Reference paper

### Yadav, R., Huser, R., and Opitz, T. (2022) A flexible Bayesian hierarchical modeling framework for spatially dependent peaks-over-threshold data. Spatial Statistics 

### This set of R codes simulate and do the inference of the product mixture model detailed below: 
Product mixture model: $Y(s)= \alpha(s) X_1(s) X_2(s) X_3(s)$,  
$X_1(s)$: i.i.d. Weibull$(1/\beta_1,1/\Gamma(1+\beta_1))$, with shape $1/\beta_1$ and scale $1/\Gamma(1+\beta_1)$, where: $\Gamma(x)$ is the gamma function  
$X_2(s)$: spatially constant term (varying with respect to time only) with marginal Weibull$(1/\beta_2,1/\Gamma(1+\beta_2))$  
$X_3(s):$ marginally Inv-Gamma$(\beta_3,\beta_3-1)$, with rate $\beta_3-1 >0$ and shape $\beta3 > 1$, dependence structure governed by Gaussian copula with exponential correlation function $\exp(-h/\rho), h\geq 0$, with range $\rho>0$  
$\alpha(s):$ scale parameter of the model with $\exp(\gamma_0 + \gamma_1 Z1(s) + \gamma_2 Z_2(s) + \gamma_3 Z_3(s))$, where $Z_i(s), i=1,2,3$ are the spatial covariates 

### The description of all the attached files are provided below:
#### Data_and_initial_values.R
R functions that simulate the data from the above product mixture model for certain combinations of model parameters and also set the initial values for all the model parameters
#### Data_nsites_20_ntime_100.Rdata
The R data that is resulting from some parameter combinations mentioned in Data_and_initial_values.R for 20 spatial locations and 100 independent temporal replicates
#### MCMC_main_function.R
This is the main MCMC function that combines several other R functions to simulate the MCMC samples of each of the model parameters
#### other-function.R
Other function that are used in MCMC_main_function.R
#### update_tuning_param_function.R
R function that tune the tuning parameters of the RWM and SGLD algorithms 
#### update.param.R
R functions that are used as a proposal distribution for all the model parameters 



