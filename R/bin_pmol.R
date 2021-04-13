#' Binned estimation of picomoles of DNA present in cfMeDIP assays
#' 
#' Given the results of model_glm_pmol and predict_pmol, adjust the predictions
#' to reflect picomoles of captured DNA overlapping a given bin in the genome. 
#' 
#' @param x     results from predict_pmol (a data.frame or GRanges) 
#'
#' @return      the same object, but with a column `adjusted_pred_con`
#'
#' @examples 
#' data(spike, package="spiky")
#' data(ssb_res, package="spiky")
#' fit <- model_glm_pmol(covg_to_df(ssb_res, spike=spike))
#' pred <- predict_pmol(fit, ssb_res, ret="df") 
#' bin_pmol(pred) 
#' 
#' @seealso model_glm_pmol
#' @seealso predict_pmol
#' 
#' @export
bin_pmol <- function(x) { 

  nameses <- names(x) 
  if (is(x, "GRanges")) nameses <- names(mcols(x))
  if (!all(c("pred_conc", "read_count") %in% nameses)) {
    stop("Please use output from predict_pmol as input to bin_pmol.")
  }
  x$adjusted_pred_con <- x$pred_conc * x$read_count
  return(x)

}
