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