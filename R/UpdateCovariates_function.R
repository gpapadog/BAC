UpdateCovariates <- function(X, Y, D, current_coefs, current_alphas,
                             current_vars, mu_priorX, mu_priorY, Sigma_priorX,
                             Sigma_priorY, omega) {
  
  num_obs <- length(X)
  num_conf <- ncol(D)
  
  # What the function will return.
  r <- list(alphas = current_alphas, coefs = current_coefs)

  for (jj in 1 : num_conf) {
    
    # ------ Updating the inclusion indicator of the exposure model. ------ #
    
    corresp_alphaY <- r$alphas[2, jj]
    
    # Prior density of coefficient = 0.
    prior_mean <- mu_priorX[jj + 1]
    prior_sd <- sqrt(diag(Sigma_priorX)[jj + 1])
    quant1 <- dnorm(x = 0, mean = prior_mean, sd = prior_sd, log = TRUE)
    
    # Prior odds of inclusion indicators.
    quant2 <- 1 / ifelse(corresp_alphaY == 0, omega + 1, 2)
    
    # Posterior density of coefficient = 0.
    rest_des_mat <- cbind(1, D[, - jj])
    resid <- X - rest_des_mat %*% current_coefs[1, - c(2, jj + 2)]
    post_quants <- PostQuants(predictor = D[, jj], resid = resid,
                              prior_mean = mu_priorX[jj + 1],
                              prior_sd = sqrt(diag(Sigma_priorX)[jj + 1]),
                              resid_var = current_vars[1])
    post_mean <- post_quants$post_mean
    post_var <- post_quants$post_var
    quant3 <- dnorm(x = 0, mean = post_mean, sd = sqrt(post_var), log = TRUE)
    
    prob1 <- exp(quant1 + quant2 - quant3)
    prob0 <- ifelse(corresp_alphaY == 0, omega / (omega + 1), 1 / 2)
    if (is.infinite(prob1) & prob1 > 0) {
      prob0 <- 0
      prob1 <- 1
    }

    r$alphas[1, jj] <- sample(c(0, 1), 1, prob = c(prob0, prob1))
    
    # ------ Updating the coefficient in the exposure model. ------ #
    
    post_sample <- rnorm(1, mean = post_mean, sd = sqrt(post_var))
    r$coefs[1, jj + 2] <- post_sample * r$alphas[1, jj]
    
    
    
    # ------ Updating the inclusion indicator of the outcome model. ------ #
    
    corresp_alphaX <- r$alphas[1, jj]

    # Prior density of coefficient = 0.
    prior_mean <- mu_priorY[jj + 2]
    prior_sd <- sqrt(diag(Sigma_priorY)[jj + 2])
    quant1 <- dnorm(x = 0, mean = prior_mean, sd = prior_sd, log = TRUE)
    
    # Prior inclusion indicator probabilities.
    quant2 <- log(ifelse(corresp_alphaX == 0, 1 / 2, omega / (omega + 1)))
    
    # Posterior density of coefficient = 0.
    rest_des_mat <- cbind(1, X, D[, - jj])
    resid <- Y - rest_des_mat %*% current_coefs[2, - (jj + 2)]
    post_quants <- PostQuants(predictor = D[, jj], resid = resid,
                              prior_mean = mu_priorY[jj + 2],
                              prior_sd = sqrt(diag(Sigma_priorY)[jj + 2]),
                              resid_var = current_vars[2])
    post_mean <- post_quants$post_mean
    post_var <- post_quants$post_var
    quant3 <- dnorm(x = 0, mean = post_mean, sd = sqrt(post_var), log = TRUE)

    prob1 <- exp(quant1 + quant2 - quant3)
    prob0 <- ifelse(corresp_alphaX == 0, 1 / 2, 1 / (omega + 1))
    if (is.infinite(prob1) & prob1 > 0) {
      prob0 <- 0
      prob1 <- 1
    }

    r$alphas[2, jj] <- sample(c(0, 1), 1, prob = c(prob0, prob1))
    
    # ------ Updating the coefficient in the outcome model. ------ #
    
    post_sample <- rnorm(1, mean = post_mean, sd = sqrt(post_var))
    r$coefs[2, jj + 2] <- post_sample * r$alphas[2, jj]
    
  }
  
  return(r)
  
}
