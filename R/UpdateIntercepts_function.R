UpdateIntercepts <- function(X, Y, D, current_coefs, current_vars, mu_priorX,
                             mu_priorY, Sigma_priorX, Sigma_priorY) {
  
  num_obs <- length(X)
  predictor <- rep(1, num_obs)
  
  # For the exposure model.
  
  resid <- X - D %*% current_coefs[1, - c(1, 2)]
  exp_post_quants <- PostQuants(predictor = predictor, resid = resid,
                                prior_mean = mu_priorX[1],
                                prior_sd = sqrt(Sigma_priorX[1, 1]),
                                resid_var = current_vars[1])
  
  resid <- Y - cbind(X, D) %*% current_coefs[2, - 1]
  out_post_quants <- PostQuants(predictor = predictor, resid = resid,
                                prior_mean = mu_priorY[1],
                                prior_sd = sqrt(Sigma_priorY[1, 1]),
                                resid_var = current_vars[2])
  
  return(c(exp_int = rnorm(1, mean = exp_post_quants$post_mean,
                           sd = sqrt(exp_post_quants$post_var)),
           out_int = rnorm(1, mean = out_post_quants$post_mean,
                           sd = sqrt(out_post_quants$post_var))))
}
