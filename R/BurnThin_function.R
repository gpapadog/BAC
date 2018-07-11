#' Burning and thinning the MCMC samples.
#' 
#' @param bac The fit from the BAC function.
#' @param burn How many of the initial samples should be dropped. Default is 0.
#' @param thin One every how many samples should be kept. Defaults to 1.
#' 
#' @example
#' bac_short <- BurnThin(bac, burn = 1000, thin = 5)
#' @export
BurnThin <- function(bac, burn = 0, thin = 1) {
  
  Nsims <- dim(bac$alphas)[3]
  if (burn >= Nsims) {
    stop('burn has to be smaller than the number of posterior samples.')
  }
  keep <- seq(burn + 1, Nsims, by = thin)
  bac_short <- list(alphas = bac$alphas[, , keep, , drop = FALSE],
                    coefs = bac$coefs[, , keep, , drop = FALSE],
                    variances = bac$variances[, , keep, drop = FALSE])
  return(bac_short)
}