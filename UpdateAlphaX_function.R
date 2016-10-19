#' Updating the outcome inclusion indicators
#' Function that updates the indicators of inclusion for predictors in the outcome
#' model.
#' 
#' @param X Vector of the treatment
#' @param D Matrix or data frame where the columns correspond to all possible
#' predictors of Y, and the rows correspond to the units.
#' @param curr_alphaX Vector of 0's and 1's of length equal to the number of columns
#' in D. Variables in D with corresponding curr_alphaX entry equal to 1 are in the
#' exposure model at the specific iteration of BAC. The proposed value of alphaX will
#' be the same as curr_alphaX by changing only one 0 to 1 or one 1 to 0.
#' @param curr_alphaY Vector of 0's and 1's of length equal to the number of columns
#' in D. Only variables in D with corresponding curr_alphaY entry equal to 1 are 
#' in the outcome model at the specific iteration of BAC.
#' @param nu_prior The value of nu in the prior of sigma square for the exposure model.
#' @param lambda_prior The value of lambda in the prior of sigma square for the
#' exposure model.
#' @param mu_prior The mean of the normal prior on the coefficients of the exposure
#' model, where the first element corresponds to the intercept, and the remaining to
#' the coefficients in from of the columns in D. The length of this vector must be
#' equal to 1 + the number of columns in D.
#' @param Sigma_prior The covariance matrix of the normal prior on the coefficients of
#' the exposure model. The structure should be in the same order as mu_prior and the
#' matrix must be symmetric with dimensions equal to 1 + the number of columns in D.
#' @param omega The omega parameter of the BAC prior.
#' 
#' @example
UpdateAlphaX <- function(X, D, curr_alphaX, curr_alphaY, nu_prior, lambda_prior,
                         mu_prior, Sigma_prior, omega) {

  prop_alphaX <- curr_alphaX
  wh_change <- sample(num_conf, 1)
  prop_alphaX[wh_change] <- 1 - prop_alphaX[wh_change]

  curr_alphaX_ind <- which(curr_alphaX == 1)
  prop_alphaX_ind <- which(prop_alphaX == 1)
  
  logAR <- CalcLogLike(outcome = X, design_mat = D[, prop_alphaX_ind],
                       nu_prior = nu_prior, lambda_prior = lambda_prior,
                       mu_prior = mu_prior[c(1, prop_alphaX_ind + 1)],
                       Sigma_prior = Sigma_prior[c(1, prop_alphaX_ind + 1),
                                                 c(1, prop_alphaX_ind + 1)]) -
    CalcLogLike(outcome = X, design_mat = D[, curr_alphaX_ind],
                nu_prior = nu_prior, lambda_prior = lambda_prior,
                mu_prior = mu_prior[c(1, curr_alphaX_ind + 1)],
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