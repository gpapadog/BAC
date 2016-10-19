#' Log Prior Odds for the BAC inclusions indicators.
#' Function that calculates the prior odds for the MH move through the model space.
#' 
#' @param proposed Vector of 0's and 1's corresponding to proposed covariates that will
#' be included in the model. 0 corresponds to exclusion and 1 to inclusion of the
#' the covariate accordingly.
#' @param current Same length with proposed. Includes the current values of the
#' inclusion indicators.
#' @param other_model Same length as proposed and current. Includes the current values
#' of the inclusion indicators for the other model. For example, when we are updating
#' the outcome model, other_model includes the inclusion indicators for the exposure.
#' @param out_model Logical. Set to TRUE if we are updating the outcome model.
#' @param omega The omega of the BAC prior.
#' 
#' @example 
LogPriorOdds <- function(proposed, current, other_model, out_model, omega) {
  wh_change <- which(proposed != current)
  r <- 1
  if (out_model) {  # If we are updating the outcome model.
    if (other_model[wh_change] == 1) {  # The exposure model includes the covariate.
      if (proposed[wh_change] == 1) {  # If the proposed model includes the covariate.
        r <- omega
      } else {  # The proposed model does not include the covariate.
        r <- 1 / omega
      }
    }
  } else {  # If we are updating the exposure model.
    if (other_model[wh_change] == 0) {
      if (proposed[wh_change] == 0) {
        r <- omega
      } else {
        r <- 1 / omega
      }
    }
  }
  return(r)
}