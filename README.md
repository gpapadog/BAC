# BAC
Bayesian Adjustment for Confounding

R package to perform the method of "Bayesian Effect Estimation Accounting for Adjustment Uncertainty" by Chi Wang, Giovanni Parmigiani, and Francesca Dominici.

# Why use BAC?

1. When the number of possible predictors of the outcome of interest is relative large,
model selection is necessary in order to choose an appropriate model.

2. The uncertainty in model selection should be propagated to effect estimates.

3. Model selection should be perfomed such that we do not miss even one important
confounder of the exposure-outcome relationship.
