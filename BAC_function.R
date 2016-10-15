BAC <- function(X, Y, D, Nsims) {
  
  num_conf <- ncol(D)
  
  alphas <- NULL
  alphas$X <- array(NA, dim = c(Nsims, num_conf))
  alphas$Y <- alphasX
  alphas$X[1, ] <- rbinom(num_conf, 1, 1 / 2)
  alphas$Y[1, ] <- rbinom(num_conf, 1, 1 / 2)
  betas <- rep(NA, Nsims)
  acc <- c(0, 0)
  
  for (ii in 2:Nsims) {
    # Update alphaX
    alphaX <- UpdateAlphaX(X = X, D = D, curr_alphaX = alphas$X[ii - 1, ],
                           curr_alphaY = alphas$Y[ii - 1, ],
                           nu_prior = nu_priorX, lambda_prior = lambda_priorX,
                           mu_prior = mu_priorX, Sigma_prior = Sigma_priorX,
                           omega = omega)
    alphas$X <- alphaX$alpha
    acc[1] <- acc[1] + alphaX$acc
    
    # Update alphaY
    alphaY <- UpdateAlphaY(Y = Y, X = X, D = D, curr_alphaX = alpha$X[ii, ],
                           curr_alphaY = alpha$Y[ii - 1, ], beta = betas[ii - 1],
                           nu_prior = nu_priorY, lambda_prior = lambda_priorY,
                           mu_prior = mu_priorY, Sigma_prior = Sigma_priorY,
                           omega = omega)
    alphas$Y <- alphaY$alpha
    acc[2] <- acc[2] + alphaY$acc
    
    # Update beta.
  }
  
  r <- NULL
  r$alphas <- alphas
  r$acc <- acc / (Nsims - 1)
  r$beta <- betas
  return(r)
}