alphas <- NULL
alphas$X <- array(NA, dim = c(Nsims, num_conf))
alphas$Y <- alphasX
alphas$X[1, ] <- rbinom(num_conf, 1, 1 / 2)
alphas$Y[1, ] <- rbinom(num_conf, 1, 1 / 2)
betas <- rep(NA, Nsims)

for (ii in 2:Nsims) {
  # Update alphaX
  curr_alphaX <- alphas$X[ii - 1, ]
  prop_alphaX <- curr_alphaX
  wh_change <- sample(num_conf, 1)
  prop_alphaX[wh_change] <- 1 - prop_alphaX[wh_change]
  curr_alphaY <- alphas$Y[ii - 1]
  
  curr_alphaX_ind <- which(curr_alphaX == 1)
  prop_alphaX_ind <- which(prop_alphaX == 1)

  logAR <- CalcLogLike(outcome = X, design_mat = D[, prop_alphaX_ind],
                       nu_prior = nu_priorX, lambda_prior = lambda_priorX,
                       mu_prior = mu_priorX[1, prop_alphaX_ind + 1],
                       Sigma_prior = Sigma_priorX[c(1, prop_alphaX_ind + 1),
                                                  c(1, prop_alphaX_ind + 1)]) -
    CalcLogLike(outcome = X, design_mat = D[, curr_alphaX_ind],
                nu_prior = nu_priorX, lambda_prior = lambda_priorX,
                mu_prior = mu_priorX[1, curr_alphaX_ind + 1],
                Sigma_prior = Sigma_priorX[c(1, curr_alphaX_ind + 1),
                                           c(1, curr_alphaX_ind + 1)])
  logAR <- logAR + logPriorOdds()
}
