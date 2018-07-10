#' MCMC update of the residual variance.
#' 
#' @param X Vector of the treatment
#' @param Y Vector of the outcome.
#' @param D Matrix or data frame where the columns correspond to all possible
#' predictors of Y, and the rows correspond to the units.
#' @param current_alphas Matrix of inclusion indicators with rows corresponding
#' to exposure and outcome models and columns to covariates.
#' @param current_coefs Matrix of coefficients with rows corresponding to
#' exposure and outcome models and columns to intercept, exposure and
#' covariates.
#' @param alpha_priorX The value of alpha in the inverse gamma prior for the
#' residual variance of the exposure model.
#' @param alpha_priorY The value of alpha in the inverse gamma prior for the
#' residual variance of the outcome model.
#' @param beta_priorX The value of beta in the inverse gamma prior for the
#' residual variance of the exposure model.
#' @param beta_priorY The value of beta in the inverse gamma prior for the
#' residual variance of the outcome model.
#' 
UpdateVariances <- function(X, Y, D, current_alphas, current_coefs,
                            alpha_priorX, alpha_priorY, beta_priorX,
                            beta_priorY) {
  
  num_obs <- nrow(X)
  r <- rep(NA, 2)
  
  # For the exposure model.
  
  alpha_new <- alpha_priorX + num_obs / 2
  des_mat <- cbind(1, D)
  resid <- X - des_mat %*% current_coefs[1, - 2]
  beta_new <- beta_priorX + sum(resid ^ 2) / 2
  r[1] <- invgamma::rinvgamma(1, shape = alpha_new, rate = beta_new)
  
  
  # For the outcome model.
  
  alpha_new <- alpha_priorY + num_obs / 2
  des_mat <- cbind(1, X, D)
  resid <- Y - des_mat %*% current_coefs[2, ]
  beta_new <- beta_priorY + sum(resid ^ 2) / 2
  r[2] <- invgamma::rinvgamma(1, shape = alpha_new, rate = beta_new)
  
  return(r)
}