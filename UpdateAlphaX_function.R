UpdateAlphaX <- function(X, D, curr_alphaX, curr_alphaY, nu_prior, lambda_prior,
                         mu_prior, Sigma_prior, omega) {

  prop_alphaX <- curr_alphaX
  wh_change <- sample(num_conf, 1)
  prop_alphaX[wh_change] <- 1 - prop_alphaX[wh_change]

  curr_alphaX_ind <- which(curr_alphaX == 1)
  prop_alphaX_ind <- which(prop_alphaX == 1)
  
  logAR <- CalcLogLike(outcome = X, design_mat = D[, prop_alphaX_ind],
                       nu_prior = nu_prior, lambda_prior = lambda_prior,
                       mu_prior = mu_prior[1, prop_alphaX_ind + 1],
                       Sigma_prior = Sigma_prior[c(1, prop_alphaX_ind + 1),
                                                 c(1, prop_alphaX_ind + 1)]) -
    CalcLogLike(outcome = X, design_mat = D[, curr_alphaX_ind],
                nu_prior = nu_prior, lambda_prior = lambda_prior,
                mu_prior = mu_prior[1, curr_alphaX_ind + 1],
                Sigma_prior = Sigma_prior[c(1, curr_alphaX_ind + 1),
                                          c(1, curr_alphaX_ind + 1)])
  logAR <- logAR + LogPriorOdds(proposed = prop_alphaX, current = curr_alphaX,
                                other_model = curr_alphaY, out_model = FALSE,
                                omega = omega)
  r <- NULL
  r$alpha <- curr_alphaX
  r$acc <- FALSE
  if (log(runif(1)) < logAR) {
    r$acc <- TRUE
    r$alpha <- prop_alphaX
  }
  return(r)
}