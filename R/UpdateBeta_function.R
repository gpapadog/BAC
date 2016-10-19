#' Updating the coefficient of treatment
#' Function that updates the coefficient in front of the treatment variable in BAC.
#' 
#' @param Y Vector of the outcome
#' @param X Vector of the treatment
#' @param D Matrix or data frame where the columns correspond to all possible
#' predictors of Y, and the rows correspond to the units.
#' @param curr_alphaY Vector of 0's and 1's of length equal to the number of columns
#' in D. Only variables in D with corresponding curr_alphaY entry equal to 1 are 
#' in the model at the specific iteration of BAC.
#' @param nu_prior The value of nu in the prior of sigma square for the outcome model.
#' @param lambda_prior The value of lambda in the prior of sigma square for the
#' outcome model.
#' @param mu_prior The mean of the normal prior on the coefficients, where the first
#' element corresponds to the intercept, the second to the treatment, and the
#' remaining to the coefficients in from of the columns in D. The length of this
#' vector must be equal to 2 + the number of columns in D.
#' @param Sigma_prior The covariance matrix of the normal prior on the coefficients.
#' The structure should be in the same order as mu_prior and the matrix must be
#' symmetric with dimensions equal to 2 + the number of columns in D.
#' 
#' @example
UpdateBeta <- function(Y, X, D, curr_alphaY, nu_prior, lambda_prior, mu_prior,
                       Sigma_prior) {
  
  n <- length(Y)
  Y <- matrix(Y, nrow = n, ncol = 1)
  curr_alphaY_ind <- which(curr_alphaY == 1)
  design_mat <- cbind(1, X, D[, curr_alphaY_ind])
  prior_Sigma <- Sigma_prior[c(1, 2, curr_alphaY_ind + 2),
                             c(1, 2, curr_alphaY_ind + 2)]
  prior_Sigma_inv <- solve(prior_Sigma)
  prior_mean <- matrix(mu_prior[c(1, 2, curr_alphaY_ind + 2)], ncol = 1,
                       nrow = length(curr_alphaY_ind) + 2)
  
  mean_betas <- solve(t(design_mat) %*% design_mat + prior_Sigma_inv) %*%
    (prior_Sigma_inv %*% prior_mean + t(design_mat) %*% Y)
  mean_betas <- matrix(mean_betas, nrow = length(mean_betas), ncol = 1)
  
  Sigma_betas <- (n + nu_prior) ^ (- 1) *
    as.numeric((nu_prior * lambda_prior +
                  t(matrix(Y - design_mat %*% mean_betas, ncol = 1, nrow = n)) %*% Y +
                  t(prior_mean - mean_betas) %*% prior_Sigma_inv %*% prior_mean)) *
    solve(t(design_mat) %*% design_mat + prior_Sigma_inv)
  
  return(rnorm(1, mean = mean_betas[2], sd = sqrt(Sigma_betas[2, 2])))
}