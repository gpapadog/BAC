#' Bayesian Adjustment for Confounding
#' Function that fits the BAC method of Wang, Parmigianni, Dominici (2012).
#' 
#' @param X Vector of the treatment
#' @param Y Vector of the outcome.
#' @param D Matrix or data frame where the columns correspond to all possible
#' predictors of Y, and the rows correspond to the units.
#' @param Nsims Number of posterior samples we wish to get.
#' 
#' @param mu_priorX The mean of the normal prior on the coefficients of the exposure
#' model, where the first element corresponds to the intercept, and the remaining to
#' the coefficients in from of the columns in D. The length of this vector must be
#' equal to 1 + the number of columns in D.
#' @param mu_priorY The mean of the normal prior on the coefficients of the outcome
#' model, where the first element corresponds to the intercept, the second to the
#' treatment, and the remaining to the coefficients in from of the columns in D. The
#' length of this vector must be equal to 2 + the number of columns in D.
#' @param nu_priorX The value of nu in the prior of sigma square for the exposure
#' model. Defaults to 0.01.
#' @param nu_priorY The value of nu in the prior of sigma square for the outcome model.
#' Defaults to 0.01.
#' @param lambda_priorX The value of lambda in the prior of sigma square for the
#' exposure model. Defaults to 0.01.
#' @param lambda_priorY The value of lambda in the prior of sigma square for the
#' outcome model. Defaults to 0.01.
#' @param Sigma_priorX The covariance matrix of the normal prior on the coefficients of
#' the exposure model. The structure should be in the same order as mu_priorX and the
#' matrix must be symmetric with dimensions equal to 1 + the number of columns in D.
#' @param Sigma_priorY The covariance matrix of the normal prior on the coefficients
#' of the outcome model. The structure should be in the same order as mu_priorY and the
#' matrix must be symmetric with dimensions equal to 2 + the number of columns in D.
#' @param omega The omega parameter of the BAC prior. Defaults to 50000.
#' 
#' @example
#' @export
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