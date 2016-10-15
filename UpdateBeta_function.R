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