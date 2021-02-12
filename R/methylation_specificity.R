#' compute methylation specificity for spike-in standards
#' 
#' In a cfMeDIP experiment, the yield of methylated fragments should be >95% 
#' (ideally 98-99%) due to the nature of the assay. 
#'
#' @param   ssb_res     output object from scan_spiked_bam
#'
#' @return              list with median and mean coverage across spike contigs
#'
#' @examples
#' data(ssb_res) 
#' methylation_specificity(ssb_res)
#' 
#' @export
methylation_specificity <- function(ssb_res) { 
  methreads <- covg_to_df(ssb_res, meth=TRUE)$read_count
  totalreads <- covg_to_df(ssb_res, meth=FALSE)$read_count
  meth_spec <- list("mean" = base::mean(methreads/totalreads), 
                    "median" = stats::median(methreads/totalreads))
  return(meth_spec)
}
