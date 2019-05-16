# BAC
Bayesian Adjustment for Confounding

R package to perform the method of "Bayesian Effect Estimation Accounting for Adjustment Uncertainty" by Chi Wang, Giovanni Parmigiani, and Francesca Dominici.

## Why use BAC?

- When the number of possible predictors of the outcome of interest is relative large,
model selection is necessary in order to choose an appropriate model.

- The uncertainty in model selection should be propagated to effect estimates.

- Model selection should be perfomed such that we do not miss even one important
confounder of the exposure-outcome relationship.

## What does BAC do?

Consider a continuous treatment X, a continuous outcome Y, and a set of possible
confounders **D** = {D<sub>1</sub>, D<sub>2</sub>, ..., D<sub>p</sub>}. BAC models:

- X given D<sub>1</sub>, D<sub>2</sub>, ..., D<sub>p</sub>

- Y given D<sub>1</sub>, D<sub>2</sub>, ..., D<sub>p</sub>

allowing only some of the D<sub>1</sub>, D<sub>2</sub>, ..., D<sub>p</sub> to be included in the exposure or the
outcome model at each iteration of the MCMC.

The BAC prior links the indicators of inclusion of D<sub>1</sub>, D<sub>2</sub>,
..., D<sub>p</sub> to encourage
the inclusion of the true confounders in the outcome model. Specifically, BAC is
designed to capture confounders that are weakly associated with the outcome, but
strongly associated with the exposure, making them strong confounders.

## What does the BAC R package do?

The BAC R package models the exposure X as normal with mean linear in the covariates
**D**, and the outcome Y also as normal with mean linear in X and in the covariates. The
priors on the coefficients are assumed normal, and the prior on the variance terms for
both the exposure and the outcome models are inverse gamma distributions.

The BAC function in the BAC R package returns a list of:
- alphaX: posterior samples of the indicators of inclusion for each of the covariates
D<sub>1</sub>, D<sub>2</sub>, ..., D<sub>p</sub> in the exposure model. Specifically
if alphaX[i, j] = 1 if in iteration i of the MCMC, Dj was included in the model,
and alphaX[i, j] = 0 otherwise.

- alphaY: posterior samples of the indicators of inclusion for each of the covariates
in the outcome model.

- beta: posterior samples of the coefficient of X in the outcome model.


## Installing BAC
Installing and using BAC in Rstudio is straightforward. You will first need the ```devtools``` R package.
### Install and load ```devtools```
Simply write ```install.packages('devtools')``` in the console to install it, and load it using ```library(devtools)```.
### Install and load ```BAC```
```
library(devtools)

devtools::install_github("gpapadog/BAC")
```


# BAC example

## toyData
toyData is a simulated data set of 100 observations with p = 50 potential confounders generated as:

- D1, D2, ..., Dp are independent N(0, 1) variables.

- X | D1, D2, ... Dp ~ N(0.1 * D1 + D2 + D3 + D4, 1)

- Y | X, D1, D2, ..., Dp ~ N(X + D1 + D2 + 0.05 * D3 + D5, 1)


## Analyzing toyData with BAC

```
library(BAC)

data(toyData)
X <- toyData[, 2]
Y <- toyData[, 1]
D <- toyData[, - c(1, 2)]

bac <- BAC(X = X, Y = Y, D = D, Nsims = 1000, chains = 3)
```
