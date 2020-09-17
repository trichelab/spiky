#' Bland-Altman plot for cfMeDIP spike standards 
#'
#' @param   fit             a model fit, from model_glm_pmol
#' 
#' @return                  a ggplot2 object
#'
#' @import  BlandAltmanLeh
#' @import  ggplot2
#' @import  scales
#' 
#' @export 
spike_bland_altman_plot <- function(fit, ...) {

  BA <- bland.altman.plot(fit$x[,1], 
                          fit$resid, 
                          conf.int = .95, 
                          col.points = fit$x[,2], 
                          pch = 19, 
                          graph.sys = "ggplot2")

  BA + 
    aes(color = fit$x[,8]) + 
    theme_bw() +
    scale_color_manual(values = c("black", "darkgrey", "navy")) +
    theme(legend.title = element_blank(), 
          legend.position = "right", 
          legend.text = element_text(size = 14), 
          axis.title.y = element_text(size = 20), 
          axis.text.x  = element_text(size = 16),
          axis.text.y = element_text(vjust = 0.5, size = 16), 
          axis.title.x = element_text(size = 20)) +
    xlab("Mean of measurements") +
    ylab("Difference")

}
