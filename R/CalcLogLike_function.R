#' Marginal log likelihood of linear model.
#' Function that calculates the marginal log likelihood of a normal outcome on a design
#' matrix with normal prior on the coefficients and inverse gamma prior on the
#' variance.
#'
#' @param outcome Vector with the entries of the outcome.
#' @param design_mat The design matrix of the outcome model.
#' @param nu_prior The value of nu in the inverse gamma prior of sigma square.
#' @param lambda_prior The value of lambda in the inverse gamma prior of sigma square.
#' @param mu_prior The mean of the normal prior on the coefficients of the model.
#' @param Sigma_prior The covariance matrix of the normal prior on the coefficients of
#' the outcome model.
CalcLogLike <- function(outcome, design_mat, nu_prior, lambda_prior, mu_prior,
                        Sigma_prior) {
  N <- length(outcome)
  if (class(design_mat) == 'numeric') {
    design_mat <- matrix(design_mat, nrow = N, ncol = 1)
  }
  if (is.data.frame(design_mat)) {
    design_mat <- as.matrix(design_mat)
  }
  design_mat <- cbind(Int = 1, design_mat)
  if (length(mu_prior) != ncol(design_mat)) {
    stop('Normal prior mean of inappropriate length.')
  }
  if (length(mu_prior) > 1) {
    if (dim(Sigma_prior)[1] != dim(Sigma_prior)[2]) {
      stop('Sigma prior has to be symmetric.')
    }
    if (dim(Sigma_prior)[1] != length(mu_prior)) {
      stop('Dimension of Sigma_prior has to be the same as the length of mu_prior.')
    }
  } 
  mu_prior <- matrix(mu_prior, ncol = 1, nrow = length(mu_prior))
  
  mat <- diag(N) + design_mat %*% Sigma_prior %*% t(design_mat)
  res_out <- outcome - design_mat %*% mu_prior
  
  r <- lgamma((nu_prior + N) / 2) + nu_prior / 2 * log(nu_prior * lambda_prior)
  r <- r - N / 2 * log(pi) - lgamma(nu_prior / 2)
  r <- r - log(det(mat) ^ (1 / 2))
  quant <- (lambda_prior * nu_prior) + t(res_out) %*% solve(mat) %*% res_out
  r <- r - (nu_prior + N) / 2 * log(quant)
  return(r)
}