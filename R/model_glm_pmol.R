#' Build a generalized linear model from spike-ins to correct bias in cfMeDIP
#'
#' formerly '2020_model_glm_fmol'.  Note that everything in x can be had from 
#' a BAM/CRAM with spike contigs named as frag_grp (len_CpGs_GC) in the index
#' FIXME: is reshape2 even used anymore?
#' 
#' @param   x     the data, with (at minimum) frag_grp, id, read_count
#' @param   conc  concentration for each spike (will be referenced if NULL)
#' @param   ...   other arguments to pass to `glm` (e.g. `family`)
#'
#' @return        the model fit for the data
#' 
#' @details       The format of 'data' is important, yet poorly documented
#'                apparently column 2 is fragment length as a character 
#'                somewhere we find the number of CpGs as an integer
#'
#' @examples 
#' 
#' data(spike_read_counts)
#' fit1 <- model_glm_pmol(spike_read_counts) 
#' 
#' data(ssb_res) # scan_spiked_bam result
#' fit2 <- model_glm_pmol(covg_to_df(ssb_res))
#' 
#' @import reshape2
#' 
#' @export
model_glm_pmol <- function(x, conc=NULL, ...) { 

  cols <- c("frag_grp", "id", "read_count")
  if (!all(cols %in% names(x))) {
    stop("Input data must have columns frag_grp, id, and read_count")
  } 

  x <- add_frag_info(x, "frag_grp") # now exported
  x$conc <- .getConcFromFraglen(x$fraglen) # added above 
  # Adjust for the 0.01ng dilution
  x$conc <- x$conc*0.9

  # Cube root CpG number to normalize
  x$CpG <- as.integer(x$CpG) # this too would be helpful to standardize 
  x$CpG_3 <- x$CpG ^ (1/3)

  # Gaussian model by default -- can alter with ... params 
  fit <- glm(formula = conc ~ read_count + fraglen + GC + CpG_3, data = x, ...)
  r2_gaussian= 1 - (fit$deviance / fit$null.deviance)
  attr(fit, "r2_gaussian") <- r2_gaussian #  report in summary
  attr(fit, "data") <- x
  return(fit) 
}

.getConcFromFraglen <- function(fraglen, concs = NULL) {
  
  if (is.null(concs)) concs <- c("80" = 0.004, "160" = 0.002, "320" = 0.001)
  res <- concs[as.character(fraglen)]
  
  # essentially the fallthrough from ifelse
  res[is.na(res)] <- concs[length(concs)]
  
  return(res)
  
}
