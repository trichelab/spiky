#' Binned estimation of picomoles of DNA present in cfMeDIP assays
#' 
#' Given the results of model_glm_pmol and predict_pmol, adjust the predictions
#' to reflect picomoles of captured DNA overlapping a given bin in the genome. 
#' 
#' @param x     results from predict_pmol (a data.frame) 
#'
#' @return      the same data.frame, but with a column `adjusted_pred_con`
#'
#' @seealso model_glm_pmol
#' @seealso predict_pmol
#' 
#' @export
bin_pmol <- function(x) {

  if (!all(c("pred_conc", "read_count") %in% names(x))) {
    stop("Please use output from predict_pmol as input")
  }
  x$adjusted_pred_con <- x$pred_conc * x$read_count
  return(x)

}
