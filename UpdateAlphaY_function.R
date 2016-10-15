UpdateAlphaY <- function(Y, X, D, curr_alphaX, curr_alphaY, beta, nu_prior,
                         lambda_prior, mu_prior, Sigma_prior, omega) {
  
  prop_alphaY <- curr_alphaY
  wh_change <- sample(num_conf, 1)
  prop_alphaY[wh_change] <- 1 - prop_alphaY[wh_change]
  
  curr_alphaY_ind <- which(curr_alphaY == 1)
  prop_alphaY_ind <- which(prop_alphaY == 1)
  
  outcome <- Y - X * beta
  
  logAR <- CalcLogLike(outcome = outcome, design_mat = D[, prop_alphaY_ind],
                       nu_prior = nu_prior, lambda_prior = lambda_prior,
                       mu_prior = mu_prior[1, prop_alphaY_ind + 1],
                       Sigma_prior = Sigma_prior[c(1, prop_alphaY_ind + 1),
                                                 c(1, prop_alphaY_ind + 1)]) -
    CalcLogLike(outcome = X, design_mat = D[, curr_alphaY_ind],
                nu_prior = nu_prior, lambda_prior = lambda_prior,
                mu_prior = mu_prior[1, curr_alphaY_ind + 1],
                Sigma_prior = Sigma_prior[c(1, curr_alphaY_ind + 1),
                                          c(1, curr_alphaY_ind + 1)])
  logAR <- logAR + LogPriorOdds(proposed = prop_alphaY, current = curr_alphaY,
                                other_model = curr_alphaX, out_model = TRUE,
                                omega = omega)
  r <- NULL
  r$alpha <- curr_alphaY
  r$acc <- FALSE
  if (log(runif(1)) < logAR) {
    r$acc <- TRUE
    r$alpha <- prop_alphaY
  }
  return(r)
}