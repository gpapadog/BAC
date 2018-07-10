#' Calculates the posterior mean and standard deviation of a coefficient.
#' 
#' @param predictor Vector including the values of the predictor for which the
#' posterior quantities are calculated.
#' @param resid The residuals of the outcome on covariates excluding intercept.
#' @param prior_mean The prior mean for the intercept.
#' @param prior_sd The prior standard deviation of the intercept.
#' @param resid_var The variance of the residuals.
#' 
PostQuants <- function(predictor, resid, prior_mean, prior_sd, resid_var) {
  
  prior_var <- prior_sd ^ 2
  
  post_var <- 1 / prior_var + sum(resid ^ 2) / resid_var
  post_var <- 1 / post_var
  
  post_mean <- prior_mean / prior_var + sum(predictor * resid) / resid_var
  post_mean <- post_var * post_mean
  
  return(list(post_mean = post_mean, post_var = post_var))
}