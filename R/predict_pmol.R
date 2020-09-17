#' predict picomoles of DNA from a fit and readcounts 
#' 
#' this seems to just wrap output from bin_pmol and model_glm_pmol... ? 
#' 
#' @param   fit result of model_glm_pmol 
#' @param   x   the data / new data 
#' 
#' @return  x, but with predictions
#' 
#' @import tools
#' @import reshape2
#' 
#' @examples 
#' 
#' data(spike_read_counts)
#' fit <- model_glm_pmol(spike_read_counts) 
#' pred <- predict_pmol(fit) 
#' 
#' @export 
predict_pmol <- function(fit, x=NULL) {

  if (is.null(x)) x <- attr(fit, "data") 
  message("start pmol prediction")
  x$pred_conc <- predict(fit, x)
  return(x) 

} 
