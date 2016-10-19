#' Updating the outcome inclusion indicators
#' Function that updates the indicators of inclusion for predictors in the outcome
#' model.
#' 
#' @param Y Vector of the outcome
#' @param X Vector of the treatment
#' @param D Matrix or data frame where the columns correspond to all possible
#' predictors of Y, and the rows correspond to the units.
#' @param curr_alphaX Vector of 0's and 1's of length equal to the number of columns
#' in D. Only variables in D with corresponding curr_alphaX entry equal to 1 are 
#' in the exposure model at the specific iteration of BAC.
#' @param curr_alphaY Vector of 0's and 1's of length equal to the number of columns
#' in D. Variables in D with corresponding curr_alphaY entry equal to 1 are in the
#' outcome model at the specific iteration of BAC. The proposed value of alphaY will be
#' the same as curr_alphaY by changing only one 0 to 1 or one 1 to 0.
#' @param beta The current value of the coefficient in front of the continuous
#' treatment.
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
#' @param omega The omega parameter of the BAC prior.
#' 
#' @example
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
                       mu_prior = mu_prior[c(1, prop_alphaY_ind + 2)],
                       Sigma_prior = Sigma_prior[c(1, prop_alphaY_ind + 2),
                                                 c(1, prop_alphaY_ind + 2)]) -
    CalcLogLike(outcome = outcome, design_mat = D[, curr_alphaY_ind],
                nu_prior = nu_prior, lambda_prior = lambda_prior,
                mu_prior = mu_prior[c(1, curr_alphaY_ind + 2)],
                Sigma_prior = Sigma_prior[c(1, curr_alphaY_ind + 2),
                                          c(1, curr_alphaY_ind + 2)])
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