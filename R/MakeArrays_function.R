#' Arrays to save MCMC posterior samples.
#' 
#' @param chains The number of MCMC chains.
#' @param Nsims The number of posterior samples per chain.
#' @param num_conf The number of potential confounders.
#' @param starting_alphas Matrix of dimensions: model (exposure or outcome),
#' potential confounders. Entries 0/1 represent exclusion/inclusion of a
#' covariate in the model. If NULL, values are set from the prior.
#' @param starting_coefs Matrix with the starting values of all coefficients.
#' Dimensions are: Exposure/Outcome model, chains, and covariate (intercept,
#' coefficient of exposure, covariates). The coefficient of exposure should be
#' NA for the exposure model. If NULL, values are set from the prior.
#' @param starting_vars Array including the starting values for the residual
#' variances. Dimensions correspond to: Exposure/Outcome model, and chains. If
#' NULL, values are set from an inverse gamma with parameters alpha and beta
#' set to the prior values times 200.
#' @param omega The omega of the BAC prior on inclusion indicators.
#' @param mu_priorX The mean of the normal prior on the coefficients of the
#' exposure model, where the first element corresponds to the intercept, and
#' the remaining to the coefficients of the columns in D.
#' @param mu_priorY The mean of the normal prior on the coefficients of the
#' outcome model, where the first element corresponds to the intercept, the
#' second to the exposure, and the remaining to the coefficients in from of the
#' columns in D.
#' @param Sigma_priorX The covariance matrix of the normal prior on the
#' coefficients of the exposure model.
#' @param Sigma_priorY The covariance matrix of the normal prior on the
#' coefficients of the outcome model.
#' @param alpha_priorX The value of alpha in the inverse gamma prior for the
#' residual variance of the exposure model.
#' @param alpha_priorY The value of alpha in the inverse gamma prior for the
#' residual variance of the outcome model.
#' @param beta_priorX The value of beta in the inverse gamma prior for the
#' residual variance of the exposure model.
#' @param beta_priorY The value of beta in the inverse gamma prior for the
#' residual variance of the outcome model.
#' 
MakeArrays <- function(chains, Nsims, num_conf, starting_alphas,
                       starting_coefs, starting_vars, omega = 50000,
                       mu_priorX, mu_priorY, Sigma_priorX, Sigma_priorY,
                       alpha_priorX, alpha_priorY, beta_priorX, beta_priorY) {
  
  
  
  # ----- Array for the inclusion indicators. ----- #
  
  alphas <- array(NA, dim = c(2, chains, Nsims, num_conf))
  dimnames(alphas) <- list(model = c('Exposure', 'Outcome'),
                           chain = 1 : chains, sample = 1 : Nsims,
                           conf = 1 : num_conf)
  
  # Starting values of the inclusion indicators.
  if (is.null(starting_alphas)) {
    starting_alphas <- array(NA, dim = c(2, chains, num_conf))
    for (cc in 1 : chains) {
      starting_alphas[2, cc, ] <- sample(c(0, 1), num_conf,
                                         prob = c(omega + 1, 2 * omega),
                                         replace = TRUE)
      for (jj in 1 : num_conf) {
        probs <- ifelse(starting_alphas[2, cc, jj] == 0, omega, 1)
        starting_alphas[1, cc, jj] <- sample(c(0, 1), 1, prob = c(probs, 1))
      }
    }
  }
  alphas[, , 1, ] <- starting_alphas
  
  
  # ----- Array for the variances. ----- #
  
  variances <- array(1, dim = c(2, chains, Nsims))
  dimnames(variances) <- list(model = c('Exposure', 'Outcome'),
                              chain = 1 : chains, sample = 1 : Nsims)
  
  # Starting values for the variances.
  if (is.null(starting_vars)) {
    starting_vars <- array(NA, dim = c(2, chains))
    starting_vars[1, ] <- invgamma::rinvgamma(chains, alpha_priorX * 200,
                                              beta_priorX * 200)
    starting_vars[2, ] <- invgamma::rinvgamma(chains, alpha_priorY * 200,
                                              beta_priorY * 200)
  }
  variances[, , 1] <- starting_vars
  
  
  # ----- Array for the coefficients. ----- #
  
  cov_names <- c('Int', 'X', paste0('C', 1 : num_conf))
  coefs <- array(0, dim = c(2, chains, Nsims, num_conf + 2))
  dimnames(coefs) <- list(model = c('Exposure', 'Outcome'),
                          chain = 1 : chains, sample = 1 : Nsims,
                          covar = cov_names)
  coefs[1, , , 2] <- NA  # No exposure coefficient for exposure model.
  
  # Starting values for the coefficients.
  if (is.null(starting_coefs)) {
    starting_coefs <- array(NA, dim = c(2, chains, num_conf + 2))
    starting_coefs[1, , - 2] <- mvnfast::rmvn(chains, mu = mu_priorX,
                                              sigma = Sigma_priorX)
    starting_coefs[2, , ] <- mvnfast::rmvn(chains, mu = mu_priorY,
                                           sigma = Sigma_priorY)
  }
  coefs[, , 1, ] <- starting_coefs
  
  return(list(alphas = alphas, coefs = coefs, variances = variances))
}