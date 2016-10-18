BAC <- function(X, Y, D, Nsims, mu_priorX = NULL, mu_priorY = NULL, nu_priorX = 0.01,
                nu_priorY = 0.01, lambda_priorX = 0.01,  lambda_priorY = 0.01,
                Sigma_priorX = NULL, Sigma_priorY = NULL, omega = 50000) {
  
  num_conf <- ncol(D)
  
  if (is.null(mu_priorX)) {
    mu_priorX <- rep(0, num_conf + 1)
  }
  if (is.null(mu_priorY)) {
    mu_priorY <- rep(0, num_conf + 2)
  }
  if (is.null(Sigma_priorX)) {
    Sigma_priorX <- 100 ^ 2 * diag(num_conf + 1)
  }
  if (is.null(Sigma_priorY)) {
    Sigma_priorY <- 100 ^ 2 * diag(num_conf + 2)
  }
  
  alphas <- NULL
  alphas$X <- array(NA, dim = c(Nsims, num_conf))
  alphas$Y <- alphas$X
  alphas$X[1, ] <- rbinom(num_conf, 1, 1 / 2)
  alphas$Y[1, ] <- rbinom(num_conf, 1, 1 / 2)
  betas <- rep(NA, Nsims)
  betas[1] <- rnorm(1, 0, 1)
  acc <- c(0, 0)
  
  for (ii in 2:Nsims) {
    
    if (ii %% 100 == 0) {
      print(paste0(round(ii / Nsims, 2) * 100, '% completed.'))
    }
    # Update alphaX
    alphaX <- UpdateAlphaX(X = X, D = D, curr_alphaX = alphas$X[ii - 1, ],
                           curr_alphaY = alphas$Y[ii - 1, ],
                           nu_prior = nu_priorX, lambda_prior = lambda_priorX,
                           mu_prior = mu_priorX, Sigma_prior = Sigma_priorX,
                           omega = omega)
    alphas$X[ii, ] <- alphaX$alpha
    acc[1] <- acc[1] + alphaX$acc
    
    # Update alphaY
    alphaY <- UpdateAlphaY(Y = Y, X = X, D = D, curr_alphaX = alphas$X[ii, ],
                           curr_alphaY = alphas$Y[ii - 1, ], beta = betas[ii - 1],
                           nu_prior = nu_priorY, lambda_prior = lambda_priorY,
                           mu_prior = mu_priorY, Sigma_prior = Sigma_priorY,
                           omega = omega)
    alphas$Y[ii, ] <- alphaY$alpha
    acc[2] <- acc[2] + alphaY$acc
    
    # Update beta.
    betas[ii] <- UpdateBeta(Y = Y, X = X, D = D, curr_alphaY = alphas$Y[ii, ],
                            nu_prior = nu_priorY, lambda_prior = lambda_priorY,
                            mu_prior = mu_priorY, Sigma_prior = Sigma_priorY)
  }
  
  r <- NULL
  r$alphas <- alphas
  r$acc <- acc / (Nsims - 1)
  r$beta <- betas
  return(r)
}