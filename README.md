# BAC
Bayesian Adjustment for Confounding

R package to perform the method of "Bayesian Effect Estimation Accounting for Adjustment Uncertainty" by Chi Wang, Giovanni Parmigiani, and Francesca Dominici.

## Why use BAC?

1. When the number of possible predictors of the outcome of interest is relative large,
model selection is necessary in order to choose an appropriate model.

2. The uncertainty in model selection should be propagated to effect estimates.

3. Model selection should be perfomed such that we do not miss even one important
confounder of the exposure-outcome relationship.

## What does BAC do?

Consider a continuous treatment X, a continuous outcome Y, and a set of possible
confounders D_1, D_2, ..., D_p. BAC models:

- X given D_1, D_2, ..., D_p

- Y given X, D_1, D_2, ..., D_p

allowing only some of the D_1, D_2, ..., D_p to be included in the exposure or the
outcome model at each iteration of the MCMC.

The BAC prior links the indicators of inclusion of D_1, D_2, ..., D_p to encourage
the inclusion of the true confounders in the outcome model. Specifically, BAC is
designed to capture confounders that are weakly associated with the outcome, but
strongly associated with the exposure, making them strong confounders.


