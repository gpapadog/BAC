UpdateExpCoef <- function(X, Y, D, current_coefs, current_vars, mu_priorY,
                          Sigma_priorY) {
  
  resid <- Y - cbind(1, D) %*% current_coefs[2, - 2]
  prior_sd <- sqrt(Sigma_priorY[2, 2])
  
  post_quants <- PostQuants(predictor = X, resid = resid,
                            prior_mean = mu_priorY[2], prior_sd = prior_sd,
                            resid_var = current_vars[2])
  
  return(rnorm(1, mean = post_quants$post_mean,
               sd = sqrt(post_quants$post_var)))
  
}