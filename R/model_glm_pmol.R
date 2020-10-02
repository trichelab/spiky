#' Attempts to build a generalized linear model from spike-in standards, 
#' to correct bias and batch effects in a cell-free MeDIP experiment. 
#'
#' formerly '2020_model_glm_fmol'.  Note that everything in X can be had from 
#' a BAM/CRAM with spike contigs named as frag_grp (len_CpGs_GC) in the index
#' 
#' @param   x     the data, with (at minimum) frag_grp, id, read_count
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
#' fit <- model_glm_pmol(spike_read_counts) 
#' 
#' @import reshape2
#' 
#' @export
model_glm_pmol <- function(x, concentration=NULL, ...) { 

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

  ##Gaussian model
  fit <- glm(formula = conc ~ read_count + fraglen + GC + CpG_3, 
             data = x, family = gaussian)
  r2_gaussian= 1 - (fit$deviance / fit$null.deviance)
  attr(fit, "r2_gaussian") <- r2_gaussian #  reprot in summary
  attr(fit, "data") <- x
  return(fit) 
}
