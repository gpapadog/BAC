#' Posterior weights of all models in MCMC.
#' 
#' Function that calculates and orders the models that were chosen in the MCMC
#' at least once.
#' 
#' @param bac The BAC fit.
#' @param model Options are 'Outcome' or 'Exposure'. Defines which model we 
#' want to look at.
#' 
#' @export
ModelWeights <- function(bac, model = c('Outcome', 'Exposure')) {
  
  num_cov <- dim(bac$alphas)[4]
  model <- match.arg(model)
  model <- ifelse(model == 'Outcome', 2, 1)
  
  alphas <- bac$alphas[model, , , ]
  alphas <- plyr::adply(alphas, 1)
  alphas <- alphas[, - 1]
  names(alphas) <- paste0('C', 1 : num_cov)
  alphas$entry <- 1 : nrow(alphas)

  alphas <- data.table::as.data.table(alphas)
  unique_alpha <- alphas[, list(number_times = length(entry)),
                         by = eval(names(alphas)[1 : num_cov])]
  unique_alpha[, proportion := number_times / nrow(alphas)]
  data.table::setorder(unique_alpha, - number_times)
  
  return(unique_alpha)
}