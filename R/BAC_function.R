#' Bayesian Adjustment for Confounding
#' 
#' Function that fits the BAC method of Wang, Parmigianni, Dominici (2012).
#' 
#' @param X Vector of the treatment
#' @param Y Vector of the outcome.
#' @param D Matrix or data frame where the columns correspond to all possible
#' predictors of Y, and the rows correspond to the units.
#' @param chains Number of MCMC chains.
#' @param Nsims Number of posterior samples we wish to get.
#' @param mu_priorX The mean of the normal prior on the coefficients of the
#' exposure model, where the first element corresponds to the intercept, and
#' the remaining to the coefficients of the columns in D. The length of this
#' vector must be equal to 1 + the number of columns in D.
#' @param mu_priorY The mean of the normal prior on the coefficients of the
#' outcome model, where the first element corresponds to the intercept, the
#' second to the exposure, and the remaining to the coefficients in from of the
#' columns in D. The length of this vector must be equal to 2 + the number of
#' columns in D.
#' @param Sigma_priorX The covariance matrix of the normal prior on the
#' coefficients of the exposure model. The dimension of the matrix should be
#' equal to 1 + the number of columns in D.
#' @param Sigma_priorY The covariance matrix of the normal prior on the
#' coefficients of the outcome model. The dimension of the matrix should be
#' equal to 2 + the number of columns in D.
#' @param alpha_priorX The value of alpha in the inverse gamma prior for the
#' residual variance of the exposure model. Defaults to 0.01.
#' @param alpha_priorY The value of alpha in the inverse gamma prior for the
#' residual variance of the outcome model. Defaults to 0.01.
#' @param beta_priorX The value of beta in the inverse gamma prior for the
#' residual variance of the exposure model. Defaults to 0.01.
#' @param beta_priorY The value of beta in the inverse gamma prior for the
#' residual variance of the outcome model. Defaults to 0.01.
#' @param omega The omega parameter of the BAC prior. Defaults to 50000.
#' @param starting_alphas Array of dimensions: model (exposure or outcome),
#' chains, potential confounders. Entries 0/1 represent exclusion/inclusion
#' of a covariate in the model. If left NULL, values are set from the prior.
#' @param starting_coefs Array with the starting values of all coefficients.
#' Dimensions are: Exposure/Outcome model, chains, and covariate (intercept,
#' coefficient of exposure, covariates). The coefficient of exposure should be
#' NA for the exposure model. If left NULL, values are set from the prior with
#' variance divided by 50 ^ 2.
#' @param starting_vars Array including the starting values for the residual
#' variances. Dimensions correspond to: Exposure/Outcome model, and chains. If
#' NULL, values are set from an inverse gamma with parameters alpha and beta
#' set to the prior values times 200.
#'  
#' @example
#' @export
BAC <- function(X, Y, D, chains, Nsims, mu_priorX = NULL, mu_priorY = NULL,
                Sigma_priorX = NULL, Sigma_priorY = NULL,
                alpha_priorX = 0.01, alpha_priorY = 0.01, beta_priorX = 0.01,
                beta_priorY = 0.01, omega = 50000, starting_alphas = NULL,
                starting_coefs = NULL, starting_vars = NULL) {
  
  D <- as.matrix(D)
  X <- matrix(X, ncol = 1)
  Y <- matrix(Y, ncol = 1)
  num_conf <- ncol(D)

  if (num_conf %in% c(0, 1)) {
    stop('If you have 0 or 1 covariate, do you need variable selection?')
  }
  
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
  
  # Checking if the dimensions of the priors are correct.
  CheckDimensions(num_conf = num_conf, mu_priorX = mu_priorX,
                  mu_priorY = mu_priorY, Sigma_priorX = Sigma_priorX,
                  Sigma_priorY = Sigma_priorY)
  
  # MCMC results' arrays, and starting values.
  arrays <- MakeArrays(chains = chains, Nsims = Nsims, num_conf = num_conf,
                       starting_alphas = starting_alphas,
                       starting_coefs = starting_coefs,
                       starting_vars = starting_vars, omega = omega,
                       mu_priorX = mu_priorX, mu_priorY = mu_priorY,
                       Sigma_priorX = Sigma_priorX,
                       Sigma_priorY = Sigma_priorY,
                       alpha_priorX = alpha_priorX,
                       alpha_priorY = alpha_priorY,
                       beta_priorX = beta_priorX, beta_priorY = beta_priorY)
  alphas <- arrays$alphas
  coefs <- arrays$coefs
  variances <- arrays$variances

  for (cc in 1 : chains) {
    for (ii in 2 : Nsims) {

      # Print progress.
      if (ii %% 100 == 0) {
        cat('Chain', cc, '-', round(ii / Nsims, 2) * 100, '% completed.\n')
      }
      
      current_alphas <- alphas[, cc, ii - 1, ]
      current_vars <- variances[, cc, ii - 1]
      current_coefs <- coefs[, cc, ii - 1, ]
      
      
      # Update the variances.
      r <- UpdateVariances(X = X, Y = Y, D = D,
                           current_alphas = current_alphas,
                           current_coefs = current_coefs,
                           alpha_priorX =  alpha_priorX,
                           alpha_priorY = alpha_priorY,
                           beta_priorX = beta_priorX,
                           beta_priorY = beta_priorY)
      variances[, cc, ii] <- r
      current_vars <- r
        
      
      # Update the inclusion indicators and coefficients.
      r <- UpdateCovariates(X = X, Y = Y, D = D, current_coefs = current_coefs,
                            current_alphas = current_alphas,
                            current_vars = current_vars, mu_priorX = mu_priorX,
                            mu_priorY = mu_priorY, Sigma_priorX = Sigma_priorX,
                            Sigma_priorY = Sigma_priorY, omega = omega)
      
      alphas[, cc, ii, ] <- r$alphas
      coefs[, cc, ii, ] <- r$coefs
      
      current_alphas <- r$alphas
      current_coefs <- r$coefs
      
      
      # Update intercepts.
      r <- UpdateIntercepts(X = X, Y = Y, D = D, current_coefs = current_coefs,
                            current_vars = current_vars, mu_priorX = mu_priorX,
                            mu_priorY = mu_priorY, Sigma_priorX = Sigma_priorX,
                            Sigma_priorY = Sigma_priorY)
      coefs[, cc, ii, 1] <- r
      current_coefs <- coefs[, cc, ii, ]
      
      
      # Update the coefficient of exposure in the outcome model.
      r <- UpdateExpCoef(X = X, Y = Y, D = D, current_coefs = current_coefs,
                         current_vars = current_vars, mu_priorY = mu_priorY,
                         Sigma_priorY = Sigma_priorY)
      coefs[2, cc, ii, 2] <- r
      current_coefs <- coefs[, cc, ii, ]
      
    }
  }
  
  return(list(alphas = alphas, coefs = coefs, variances = variances))
  
}
