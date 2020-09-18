#' Bland-Altman plot for cfMeDIP spike standards 
#'
#' @param   fit             a model fit, from predict_pmol (?)
#' 
#' @return                  a ggplot2 object
#'
#' @import  BlandAltmanLeh
#' @import  ggplot2
#' @import  scales
#'
#' @examples 
#' 
#' data(spike_read_counts)
#' fit <- model_glm_pmol(spike_read_counts) 
#' spike_bland_altman_plot(fit) 
#' 
#' @export 
spike_bland_altman_plot <- function(fit) { 

  if (!is.data.frame(fit$x)) fit$x <- attr(fit, "data") 
  if (!"pred_conc" %in% names(fit$x)) fit$x <- predict_pmol(fit) 
  fit$x$fragment_len <- as.factor(fit$x$fragment_len)
  
  BA <- bland.altman.plot(fit$x$pred_conc,
                          fit$resid, 
                          conf.int = .95, 
                          col.points = fit$x[,2], 
                          pch = 19, 
                          graph.sys = "ggplot2")

  # what is fit$x[,8] ? ##This is the fragment length variable
  BA + 
    aes(color = fit$x$fragment_len) +  # ? 
    theme_bw() +
    scale_color_manual(values = c("black", "darkgrey", "navy")) +
    theme(legend.title = element_blank(), 
          legend.position = "right", 
          legend.text = element_text(size = 14), 
          axis.title.y = element_text(size = 20), 
          axis.text.x  = element_text(size = 16),
          axis.text.y = element_text(vjust = 0.5, size = 16), 
          axis.title.x = element_text(size = 20)) +
    xlab("Mean of measurements (picomoles)") +
    ylab("Difference (picomoles)")

}
