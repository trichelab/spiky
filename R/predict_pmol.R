#' predict picomoles of DNA from a fit and readcounts 
#' 
#' this seems to just wrap output from bin_pmol and model_glm_pmol... ? 
#' 
#' @param   fit result of model_glm_pmol 
#' @param   x   the data
#' 
#' @return  predictions
#' 
#' @import tools
#' @import reshape2
#' 
#' @export 
predict_pmol <- function(fit, x=NULL) {

  if (is.null(x)) x <- attr(fit, "data") 

  # probably unnecessary, revisit later 
  data_melt <- melt(x,
                    id.vars = c("GC", "UMI", "fraglen", "CpG", "pos", "chr"),
                    measure.vars = "read_count")

  # all of this is probably unnecessary if x is extracted from fit
  data_melt$value <- as.numeric(as.character(data_melt$value))
  data_melt$CpG_3 <- (data_melt$CpG) ^ (1 / 3) # should already be in place? 
  # typically a fit object carries around its own data in a slot
  # if not, we can just attr() the data into it
  data_melt$GC <- as.numeric(as.character(data_melt$GC))
  data_melt$fragment_len <- as.numeric(data_melt$fraglen)


  print("start fmol prediction")
  ##Predict fmol values
  pred_conc <- as.data.frame(predict(fit, data_melt))
  names(pred_conc) <- "value"
  ##Replace value with conc
  data_melt$value <- NULL
  data_melt <- cbind(data_melt, pred_conc)
  data_melt$variable <- NULL
  colnames(data_melt) <- c("GC","UMI", "fragment_len", "CpG", "pos", "chr", "CpG3", "fmol")
  data_melt <- data_melt[,c("chr", "pos", "fragment_len", "UMI", "CpG", "CpG3", "GC", "fmol")]

  # print(head(data_melt))
  return(data_melt) 

} 
