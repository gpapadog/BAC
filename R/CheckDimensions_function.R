CheckDimensions <- function(num_conf, mu_priorX, mu_priorY, Sigma_priorX,
                            Sigma_priorY) {
  
  if (any(length(mu_priorX) != num_conf + 1 |
          length(mu_priorY) != num_conf + 2 |
          nrow(Sigma_priorX) != ncol(Sigma_priorX) | 
          nrow(Sigma_priorY) != ncol(Sigma_priorY) |
          nrow(Sigma_priorX) != num_conf + 1 | 
          nrow(Sigma_priorY) != num_conf + 2)) {
    
    stop('Priors on coefficients are not of correct dimensions.')
    
  }
  
}