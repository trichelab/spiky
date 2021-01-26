#' compute methylation specificity for spike-in standards
#' 
#' In a cfMeDIP experiment, the yield of methylated fragments should be >95% 
#' (ideally 98-99%) due to the nature of the assay. 
#'
#' @export
methylation_specificity <- function(ssb_res, ...) { 
  methreads <- covg_to_df(ssb_res, meth=TRUE)$read_count
  totalreads <- covg_to_df(ssb_res, meth=FALSE)$read_count
  meth_spec <- list("mean" = mean(methreads/totalreads), 
                    "median" = median(methreads/totalreads))
  return(meth_spec)
  # data("spike")
  # methylated_spikes <- rownames(spike[spike$methylated == 1,])
  # sum(spike_data[methylated_spikes])
  # methyl_spec <- sum(spike_data[methylated_spikes]) / (sum(spike_data))
  # return(methyl_spec)
}
